## GB_READ: READ A GENBANK FILE IN RUST

This code uses the [`gb-io`](https://github.com/dlesl/gb-io) crate for parsing a GenBank-formatted file. 

To use:
1) Clone this GitHub repository (`git clone https://github.com/jperkel/gb_read.git`).
2) Build in Rust with `cargo run`. (If you do not have gb-io installed, this step will do so for you.)
3) By default, the program will parse the included file `nc_005816.gb`. However, the user can specify another filename at the command line, eg: `cargo run pbr322.gb`.

