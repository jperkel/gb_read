#[macro_use]
extern crate gb_io;
extern crate regex;

use std::cmp;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io;
use std::io::Write;
use std::process;

use gb_io::reader::SeqReader;
use regex::Regex;

// Translation table 11 from https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi.
// Indexing into this array is as follows:
// T = 0, C = 1, A = 2, G = 3
// index = (16 * first_base) + (4 * second_base) + third_base
// So... TTT = 0, TTC = 1, TTA = 2, ... , GGC = 61, GGA = 62, GGG = 63
const GENETIC_CODE: &str = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

// ERR_BAD_NT is an error value for an invalid nucleotide
const ERR_BAD_NT: usize = 99;


// map a base for indexing the GENETIC_CODE string
// x (char): base to look up
// returns: usize
fn lookup(x: char) -> usize {
    match x {
        'T' => 0,
        'C' => 1,
        'A' => 2,
        'G' => 3,
        _ => ERR_BAD_NT, // unknown base
    }
}

// get the three-letter equivalent of an amino acid
// aa (char): the amino acid to translate, eg 'A' -> "Ala"
// returns: String
fn three_letter_code(aa: char) -> String {
    let three_letter_map: HashMap<char, &str> = [
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
    .collect();

    // check input.
    let re = Regex::new(r"[A-Z*]").unwrap();
    if !re.is_match(&aa.to_string()) {
        println!("Invalid amino acid: {}", aa);
        process::exit(1);
    }
    three_letter_map[&aa].to_string()
}

// translate a codon into its corresponding amino acid
// triplet (&str): a three-letter codon eg "ATG"
// returns: char
fn translate(triplet: &str) -> char {
    let mut codon = vec![ERR_BAD_NT; 3];

    for (i, base) in triplet.chars().enumerate() {
        codon[i] = lookup(base);
    }

    if codon.contains(&ERR_BAD_NT) {
        println!("Something went wrong in translation!");
        println!("{:?}", codon);
        process::exit(1);
    }

    let index: usize = (codon[0] * 16) + (codon[1] * 4) + codon[2];
    // translate the codon into single-letter code
    let c = GENETIC_CODE.chars().nth(index).unwrap();

    c
}

// print a pretty DNA sequence and its translation, plus line numbering
// eg: 001 MetSerIle...
//     001 ATGAGTATT... 
//
// s (&str): DNA sequence to print
fn print_seq(s: &str) {
    let linelen = 72; // print 72 bases per line (24 amino acids)

    // how many lines to print
    let mut nlines = s.len() / linelen;
    // if there's a remainder, add one line.
    if s.len() % linelen != 0 {
        nlines += 1;
    }

    // how wide does the numbering block need to be?
    let ndigits = count_digits(s.len() as u16);

    // get the translation of this sequence
    let mut peptide3 = String::new();
    let n_codons = s.len() / 3;
    for i in 0..n_codons {
        let codon = &s[i * 3..(i * 3) + 3]; // take a 3-base slice of the sequence
        let aa = translate(&codon);
        // translate and add to the string
        peptide3.push_str(&three_letter_code(aa));
    }

    for i in 0..nlines {
        let begin = i * linelen;
        // adjust 'end' if near the end of the sequence
        let end = cmp::min((i * linelen) + linelen, s.len());
        
        // print translation
        println!(
            "{number:>0width$} {}",
            &peptide3[begin..end],
            number = (begin / 3) + 1, // divide by 3 b/c 3 bases/amino acid
            width = ndigits
        );

        // print DNA
        println!(
            "{number:>0width$} {}\n",
            &s[begin..end],
            number = (begin) + 1,
            width = ndigits
        );
    }
}

// count the digits in a number
// n (u16): number to count
// returns: usize
//
// NOTE: n is type u16, so allowable input is 0..65535.
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
    let args: Vec<String> = env::args().collect();
    if args.len() > 1 {
        filename = &args[1];
    }
    if !std::path::Path::new(filename).exists() {
        println!("File '{}' does not exist.", filename);
        process::exit(1);
    }

    println!("\nReading records from file '{}'...", filename);

    let file = File::open(filename).unwrap();

    // genes, descs and locs are vectors of Strings to hold gene names,
    // gene descriptions, and gene locations, respectively
    let mut genes = Vec::<String>::new();
    let mut descs = Vec::<String>::new();
    let mut locs = Vec::<gb_io::seq::Location>::new();

    // a HashMap holding counts of each feature type, eg, "("CDS", 5)"
    let mut feature_map: HashMap<String, usize> = HashMap::new();

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

                locs.push(f.location.clone());
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
            let mut retval = String::new();
            io::stdin()
                .read_line(&mut retval)
                .expect("Failed to read from stdin");

            // use trim() to delete the trailing newline ('\n') char
            retval = retval.trim().to_string();
            if (retval == "q") || (retval == "Q") {
                process::exit(0);
            }

            let selection = match retval.parse::<usize>() {
                Ok(i) => i, // if good input, just return the number
                Err(_) => {
                    println!("Invalid input: '{}'", retval);
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
                    let s =
                        String::from_utf8(seq.extract_location(&locs[selection]).unwrap().to_vec())
                            .unwrap()
                            .to_ascii_uppercase();

                    println!("\n{}: {}", genes[selection], descs[selection]);
                    print_seq(&s);
                    println!("DNA:     {:>5} bases", s.len());
                    println!("Protein: {:>5} amino acids\n", s.len()/3)
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
        assert_eq!(translate("ATG"), 'M');
    }

    #[test]
    fn test_translate_tag() {
        assert_eq!(translate("TAG"), '*');
    }

    #[test]
    fn test_translate_ttt() {
        assert_eq!(translate("TTT"), 'F');
    }

    #[test]
    fn test_one_to_three_translate() {
        assert_eq!(three_letter_code(translate("ATG")), "Met");
    }

    #[test]
    fn test_translate_pyr() {
        assert_eq!(three_letter_code('O'), "Pyr");
    }

    #[test]
    fn test_translate_sel() {
        assert_eq!(three_letter_code('U'), "Sel");
    }

    #[test]
    fn translate_bad_aa() {
        assert_eq!(three_letter_code('J'), "???");
    }

    #[test]
    fn test_count_digits() {
        assert_eq!(count_digits(10000), 5);
    }

    #[test]
    fn test_count_digits2() {
        assert_eq!(count_digits(2500), 4);
    }
}
