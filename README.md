## GB_READ: READ A GENBANK FILE IN RUST

This code uses the [`gb-io`](https://github.com/dlesl/gb-io) crate for parsing a GenBank-formatted file, counting the genes, and translating them. 

To use:
1) If you have not installed the Rust programming language, do so [here](https://www.rust-lang.org/tools/install). 
2) Clone this GitHub repository (`git clone https://github.com/jperkel/gb_read.git`).
3) Build the application with `cargo run`. (If you do not have gb-io installed, this step will do so for you.)
4) By default, the program will parse the included file `nc_005816.gb`. However, the user can specify another filename at the command line, eg: `cargo run pbr322.gb`.
