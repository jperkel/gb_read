use std::cmp;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io;
use std::io::Write;
use std::process;

use gb_io::reader::SeqReader;
use gb_io::{feature_kind, qualifier_key};
use once_cell::sync::Lazy;

#[derive(thiserror::Error, Debug)]
enum Error {
    #[error("Something went wrong in translation! {0:?}")]
    BadNucleotide(Vec<usize>),

    #[error("Invalid amino acid: {0}")]
    InvalidAminoAcid(char),
}

/// Translation table 11 from https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi.
/// Indexing into this array is as follows:
/// T = 0, C = 1, A = 2, G = 3
/// index = (16 * first_base) + (4 * second_base) + third_base
/// So... TTT = 0, TTC = 1, TTA = 2, ... , GGC = 61, GGA = 62, GGG = 63
const GENETIC_CODE: &[u8] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

/// some codons can be used as Met in the start position
const STARTS: &[u8] = b"---M------**--*----M------------MMMM---------------M------------";

/// ERR_BAD_NT is an error value for an invalid nucleotide
const ERR_BAD_NT: usize = 99;

/// map a base for indexing the GENETIC_CODE string
/// x (char): base to look up
/// returns: usize
const fn lookup(x: u8) -> usize {
    match x as char {
        'T' => 0,
        'C' => 1,
        'A' => 2,
        'G' => 3,
        _ => ERR_BAD_NT, // unknown base
    }
}

static THREE_LETTER_CODE: Lazy<HashMap<char, &'static str>> = Lazy::new(|| {
    [
        ('A', "Ala"),
        ('B', "???"),
        ('C', "Cys"),
        ('D', "Asp"),
        ('E', "Glu"),
        ('F', "Phe"),
        ('G', "Gly"),
        ('H', "His"),
        ('I', "Ile"),
        ('J', "???"),
        ('K', "Lys"),
        ('L', "Leu"),
        ('M', "Met"),
        ('N', "Asn"),
        ('O', "Pyr"),
        ('P', "Pro"),
        ('Q', "Gln"),
        ('R', "Arg"),
        ('S', "Ser"),
        ('T', "Thr"),
        ('U', "Sel"),
        ('V', "Val"),
        ('W', "Trp"),
        ('X', "???"),
        ('Y', "Tyr"),
        ('Z', "???"),
        ('*', "***"),
    ]
    .iter()
    .copied()
    .collect()
});

/// get the three-letter equivalent of an amino acid
/// aa (char): the amino acid to translate, eg 'A' -> "Ala"
/// returns: String
fn three_letter_code(aa: char) -> Result<String, Error> {
    // check input.
    match aa {
        'A'..='Z' => Ok(THREE_LETTER_CODE[&aa].to_string()),
        _ => Err(Error::InvalidAminoAcid(aa)),
    }
}

/// translate a codon into its corresponding amino acid
/// triplet (&str): a three-letter codon eg "ATG"
/// i (usize): codon position. if 0, use the STARTS table
/// returns: char
fn translate(triplet: &[u8], i: usize) -> Result<char, Error> {
    let mut codon = vec![ERR_BAD_NT; 3];

    for (i, base) in triplet.iter().enumerate() {
        codon[i] = lookup(*base);
    }

    if codon.contains(&ERR_BAD_NT) {
        return Err(Error::BadNucleotide(codon));
    }

    let index: usize = (codon[0] * 16) + (codon[1] * 4) + codon[2];
    // translate the codon into single-letter code

    let c = if (i == 0) && (STARTS[index] == b'M') {
        b'M'
    } else {
        GENETIC_CODE[index]
    };

    Ok(c as char)
}

/// print a pretty DNA sequence and its translation, plus line numbering
/// eg: 001 MetSerIle...
///     001 ATGAGTATT...
///
/// s (&str): DNA sequence to print
fn print_seq(s: &str) -> Result<(), Error> {
    let line_len = 72; // print 72 bases per line (24 amino acids)

    // how many lines to print
    let n_lines = if s.len() % line_len != 0 {
        (s.len() / line_len) + 1
    } else {
        s.len() / line_len
    };

    // how wide does the numbering block need to be?
    let n_digits = count_digits(s.len() as u16);

    // get the translation of this sequence
    let mut peptide = String::new();
    for (i, codon) in s.as_bytes().windows(3).enumerate() {
        let aa = translate(codon, i)?;
        // translate and add to the string
        peptide.push_str(&three_letter_code(aa)?);
    }

    for i in 0..n_lines {
        let begin = i * line_len;
        // adjust 'end' if near the end of the sequence
        let end = cmp::min((i * line_len) + line_len, s.len());

        // print translation
        println!(
            "{number:>0width$} {}",
            &peptide[begin..end],
            number = (begin / 3) + 1, // divide by 3 b/c 3 bases/amino acid
            width = n_digits
        );

        // print DNA
        println!(
            "{number:>0width$} {}\n",
            &s[begin..end],
            number = (begin) + 1,
            width = n_digits
        );
    }

    Ok(())
}

/// count the digits in a number
/// n (u16): number to count
/// returns: usize
///
/// NOTE: n is type u16, so allowable input is 0..65535.
fn count_digits(mut n: u16) -> usize {
    let mut digits: usize = 1;

    while n > 0 {
        n /= 10;
        if n >= 1 {
            digits += 1;
        }
    }

    digits
}

fn main() {
    // file to process; user can override on the command line
    let mut filename = "nc_005816.gb";

    // check to see if user provided an alternate file name...
    let args: Vec<_> = env::args().collect();
    if args.len() > 1 {
        filename = &args[1];
    }
    if !std::path::Path::new(filename).exists() {
        println!("File '{}' does not exist.", filename);
        process::exit(1);
    }

    println!("\nReading records from file '{}'...", filename);

    let file = File::open(filename).unwrap();

    // vectors to hold gene names and gene descriptions
    let mut genes = vec![];
    let mut descs = vec![];

    // HashMap holding counts of each feature type, eg, "("CDS", 5)"
    let mut feature_map = HashMap::new();

    // length of the longest feature type, for printing
    let mut feat_len = 0;

    for r in SeqReader::new(file) {
        let seq = r.unwrap();
        let mut gene_count = 0;
        let mut feat_count = 0;
        let record_name = seq.name.clone().unwrap();
        println!("Record name: {}", record_name);
        println!("Sequence length: {}", seq.len());

        for f in &seq.features {
            feat_count += 1;

            // count the different feature types. if a type hasn't been seen
            // before, set its value to zero and add 1
            *feature_map.entry(f.kind.to_string()).or_insert(0) += 1;

            if f.kind.to_string().len() > feat_len {
                feat_len = f.kind.to_string().len();
            }

            // collect protein_id, product, and location data for each "CDS", ie gene
            if f.kind == feature_kind!("CDS") {
                let gene = f
                    .qualifier_values(qualifier_key!("protein_id"))
                    .next()
                    .unwrap()
                    .to_string()
                    .replace('\n', "");
                genes.push(gene);

                let desc = f
                    .qualifier_values(qualifier_key!("product"))
                    .next()
                    .unwrap()
                    .to_string()
                    .replace('\n', "");
                descs.push(desc);

                gene_count += 1;
            }
        }

        println!(
            "\nFound {} features, including {} genes.",
            feat_count, gene_count
        );

        for (key, count) in feature_map.iter() {
            println!("{k:<0l$}: {}", count, k = key, l = feat_len + 1);
        }

        println!();
        for i in 0..gene_count {
            println!("{}) {}: {}", i, genes[i], descs[i]);
        }

        if gene_count > 0 {
            let range = {
                if gene_count == 1 {
                    "0".to_string()
                } else {
                    format!("0-{}", gene_count - 1)
                }
            };
            print!("\nSelect a gene number [{}] or 'q' to quit: ", range);
            // flush buffer...
            io::stdout().flush().unwrap();

            // get input...
            let mut ret_val = String::new();
            io::stdin()
                .read_line(&mut ret_val)
                .expect("Failed to read from stdin");

            // use trim() to delete the trailing newline ('\n') char
            ret_val = ret_val.trim().to_string();
            if (ret_val == "q") || (ret_val == "Q") {
                process::exit(0);
            }

            let selection = match ret_val.parse::<usize>() {
                Ok(i) => i, // if good input, just return the number
                Err(_) => {
                    println!("Invalid input: '{}'", ret_val);
                    process::exit(1);
                }
            };

            if selection > gene_count - 1 {
                println!("Invalid input: '{}'", selection);
                process::exit(1);
            }

            // find the requested gene and print it.
            for f in &seq.features {
                if f.kind == feature_kind!("CDS")
                    && f.qualifier_values(qualifier_key!("protein_id"))
                        .next()
                        .unwrap()
                        == genes[selection]
                {
                    let s = String::from_utf8(
                        seq.extract_location(&f.location.clone()).unwrap().to_vec(),
                    )
                    .unwrap()
                    .to_ascii_uppercase();

                    println!("\n{}: {}", genes[selection], descs[selection]);
                    print_seq(&s).unwrap();
                    println!("DNA:     {:>5} bases", s.len());
                    println!("Protein: {:>5} amino acids\n", s.len() / 3)
                }
            }
        }
    }
}

// tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_translate_atg() {
        assert!(matches!(translate(b"ATG", 1), Ok('M')));
    }

    #[test]
    fn test_translate_atg_as_start() {
        assert!(matches!(translate(b"ATG", 0), Ok('M')));
    }

    #[test]
    fn test_translate_gtg() {
        assert!(matches!(translate(b"GTG", 1), Ok('V')));
    }

    #[test]
    fn test_translate_gtg_as_start() {
        assert!(matches!(translate(b"GTG", 0), Ok('M')));
    }

    #[test]
    fn test_translate_tag() {
        assert!(matches!(translate(b"TAG", 1), Ok('*')));
    }

    #[test]
    fn test_translate_ttt() {
        assert!(matches!(translate(b"TTT", 1), Ok('F')));
    }

    #[test]
    fn test_one_to_three_translate() {
        assert_eq!(
            three_letter_code(translate(b"ATG", 0).unwrap()).unwrap(),
            "Met"
        );
    }

    #[test]
    fn test_translate_pyr() {
        assert_eq!(three_letter_code('O').unwrap(), "Pyr");
        //assert!(matches!(three_letter_code('O'), Ok("Pyr")));
    }

    #[test]
    fn test_translate_sel() {
        assert_eq!(three_letter_code('U').unwrap(), "Sel");
    }

    #[test]
    fn translate_bad_aa() {
        assert_eq!(three_letter_code('J').unwrap(), "???");
    }

    #[test]
    fn test_count_digits1() {
        assert_eq!(count_digits(10000), 5);
    }

    #[test]
    fn test_count_digits2() {
        assert_eq!(count_digits(2500), 4);
    }
}
