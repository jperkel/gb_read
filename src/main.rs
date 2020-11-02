#[macro_use]
extern crate gb_io;

use std::fs::File;
use std::collections::HashMap;
use std::io;
use std::io::Write;
use std::process;
use std::env;

use gb_io::reader::SeqReader;

// Translation table 11 from https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. 
// Indexing into this array is as follows: 
// T = 0, C = 1, A = 2, G = 3
// index = (16 * first_base) + (4 * second_base) + third_base 
// So... TTT = 0, TTC = 1, TTA = 2, ... , GGC = 61, GGA = 62, GGG = 63
const GENETIC_CODE: &str = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

// ERR_BAD_NT is an error value for an invalid nucleotide 
const ERR_BAD_NT: usize = 99; 

// enumeration of genetic sequence types
enum SeqType {
    DNA,
    Protein1, // single-letter amino acid code
    Protein3, // triple-letter amino acid code
}

// enumeration of translation types
enum Translation {
    OneLetter,
    ThreeLetter,
}

// given an input 'char', return a base equivalent
fn lookup(x: char) -> usize {
    match x {
        'T' => return 0,
        'C' => return 1,
        'A' => return 2,
        'G' => return 3,
        _ => return ERR_BAD_NT, // unknown base
    }
}

// translate a codon into its corresponding amino acid
fn translate(triplet: &str, t: Translation) -> String {
    let three_letter_code: HashMap<char, &str> = [
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
    ].iter().copied().collect();

    let mut codon = vec![ERR_BAD_NT; 3];

    for (i,base) in triplet.chars().enumerate() {
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
    match t {
        Translation::OneLetter => return c.to_string(),
        Translation::ThreeLetter => return three_letter_code[&c].to_string(),
    }
}

// print a pretty sequence, 72 bases per line, plus base numbering
// s: sequence
// t: sequence type (DNA, Protein1 or Protein3)
fn print_seq(s: &str, t: SeqType) {
    let linelen = 72;

    let divisor = match t {
        // if we're printing a Protein3, count amino acids, not bases, 
        // so divide by 3
        SeqType::DNA | SeqType::Protein1 =>  1,
        SeqType::Protein3 => 3,
    };

    // how many lines to print
    let mut nlines = s.len()/linelen;
    // if there's a remainder, add one line.
    if s.len() % linelen != 0 {
        nlines += 1;
    }

    // print the lines
    for i in 0..nlines {
        let start = i*linelen;
        let mut end = (i*linelen) + linelen;
        if end > s.len() {
            end = s.len();
        }
        let myline = &s[start..end];
        println!("{number:>0width$} {}", myline, number=((start/divisor))+1, width=count_digits(s.len() as u16));
    }
}

// NOTE: n is type u16, so allowable input is 0..65535. 
fn count_digits(mut n: u16) -> usize {
    let mut digits: usize = 1;

    while n > 0 {
        n = n/10;
        if n >= 1 {
            digits += 1;
        }
    };

    return digits;
}

fn main() {
    let mut filename = "nc_005816.gb";

    // the user can provide another file on the command line
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
    let mut genes = Vec::<String>::new();
    let mut descs = Vec::<String>::new();
    let mut locs = Vec::<gb_io::seq::Location>::new();

    for r in SeqReader::new(file) {
        let seq = r.unwrap();
        let mut gene_count = 0;
        let mut feat_count = 0;
        let record_name = seq.name.clone().unwrap();
        println!("Record name: {}", record_name);
        println!("Sequence length: {}", seq.len());

        for f in &seq.features {
            feat_count += 1;
            // collect protein_id, product, and location data for each "CDS", ie gene
            if f.kind==feature_kind!("CDS") {
                let gene = f.qualifier_values(qualifier_key!("protein_id"))
                    .next()
                    .unwrap()
                    .to_string()
                    .replace('\n',"");
                genes.push(gene);

                let desc = f.qualifier_values(qualifier_key!("product"))
                    .next()
                    .unwrap()
                    .to_string()
                    .replace('\n',"");
                descs.push(desc);

                locs.push(f.location.clone());
                gene_count += 1;
            }
        }
        println!("\nFound {} features, including {} genes.", feat_count, gene_count);
        for i in 0..gene_count {
            println!("{}) {}: {}", i, genes[i],descs[i]);
        }

        if gene_count > 0 {
            print!("\nWhich would you like to view [0-{}]: ", gene_count-1);
            // flush buffer...
            io::stdout().flush().unwrap();

            // get input...
            let mut retval = String::new();
            io::stdin().read_line(&mut retval).expect("Failed to read from stdin");

            // use trim() to delete the trailing newline ('\n') char
            retval = retval.trim().to_string();

            let selection = match retval.parse::<usize>() {
                Ok(i) => i, // if good input, just return the number
                Err(_) => {
                    println!("Invalid input: '{}'", retval);
                    process::exit(1);
                },
            };

            if selection > gene_count-1 {
                println!("Invalid input: '{}'", selection);
                process::exit(1);
            }

            println!("You selected: {}", selection);

            for f in &seq.features {
                if f.kind==feature_kind!("CDS") {
                    if f.qualifier_values(qualifier_key!("protein_id")).next().unwrap() == genes[selection] {
                        let s = String::from_utf8(seq.extract_location(&locs[selection]).unwrap().to_vec())
                            .unwrap()
                            .to_ascii_uppercase();

                        println!("\nDNA sequence:");
                        print_seq(&s, SeqType::DNA);
                        println!("Length: {}\n", s.len());
        
                        let mut peptide1 = String::new();    
                        let mut peptide3 = String::new();    
                        let n_codons = s.len()/3;
                        for i in 0..n_codons {
                            let codon = &s[i*3..(i*3)+3]; // take a 3-base slice of the sequence
                            peptide1.push_str(&translate(&codon,    Translation::OneLetter)); // translate and add to the string
                            peptide3.push_str(&translate(&codon, Translation::ThreeLetter)); // translate and add to the string
                        }
                        println!("One-letter code:");
                        print_seq(&peptide1, SeqType::Protein1);
                        println!("\nThree-letter code:");
                        print_seq(&peptide3, SeqType::Protein3);
                        println!("Length: {} (including stop)\n", n_codons);
                    }
                }
            }
        }
    }
}