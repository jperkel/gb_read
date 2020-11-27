## GB_READ: READ A GENBANK FILE IN RUST

Update 26 Nov 2020: Thanks to a suggestion from [Avi Srivastava](https://github.com/k3yavi), gb_read now uses the [`clap`](https://crates.io/crates/clap) crate for command-line processing. As a result, the syntax for specifying an alternate file on the command line has changed. See step 4 below for details.

This code uses the [`gb-io`](https://github.com/dlesl/gb-io) crate for parsing a GenBank-formatted file, counting the genes, and translating them. 

To use:
1) If you have not installed the Rust programming language, do so [here](https://www.rust-lang.org/tools/install). 
2) [Clone this GitHub repository](https://docs.github.com/en/free-pro-team@latest/github/creating-cloning-and-archiving-repositories/cloning-a-repository). Click the green "Code" button above the file listing and select "Open with GitHub Desktop". Or, from the command line, execute `git clone https://github.com/jperkel/gb_read.git`.
3) Build the application. From the command line, execute `cargo run`. (If you do not have gb-io installed, this step will do so for you.)
4) By default, the program will parse the included file [`nc_005816.gb`](https://github.com/jperkel/gb_read/blob/main/nc_005816.gb). However, the user can specify another filename at the command line, eg: ~~`cargo run pbr322`~~[`cargo run -- -i pbr322.gb`](https://github.com/jperkel/gb_read/blob/main/pbr322.gb).
