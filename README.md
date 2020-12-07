# GB_READ: A Rust GenBank file reader

**Update: 26 Nov 2020**: Thanks to a suggestion from [Avi Srivastava](https://github.com/k3yavi), gb_read now uses the [`clap`](https://crates.io/crates/clap) crate for command-line processing. As a result, the syntax for specifying an alternate file on the command line has changed. See step 4 below for details.

This project was created to accompany a [*Nature* Toolbox article](https://www.nature.com/articles/d41586-020-03382-2) (published 1 Dec 2020) on the use of Rust in science. It uses the [`gb-io`](https://github.com/dlesl/gb-io) crate for parsing a GenBank-formatted file, counting the genes, and translating them. 

To use:
1) If you have not installed the Rust programming language, do so [here](https://www.rust-lang.org/tools/install). 
2) [Clone this GitHub repository](https://docs.github.com/en/free-pro-team@latest/github/creating-cloning-and-archiving-repositories/cloning-a-repository). Click the green "Code" button above the file listing and select "Open with GitHub Desktop". Or, from the command line, execute `git clone https://github.com/jperkel/gb_read.git`.
3) Build the application, using Rust's build tool/package manager, Cargo. From the command line, execute `cargo run`. (If you do not have gb-io and other required libraries (or 'crates') installed, this step will do so for you.)
4) By default, the program will parse the included file [`nc_005816.gb`](https://github.com/jperkel/gb_read/blob/main/nc_005816.gb). However, the user can specify another filename at the command line, eg: ~~`cargo run pbr322`~~`cargo run -- -i pbr322.gb`. ([`pbr322.gb`](https://github.com/jperkel/gb_read/blob/main/pbr322.gb) and two additional example files, [`puc19.gb`](https://github.com/jperkel/gb_read/blob/main/puc19.gb) and [`circdna.gb`](https://github.com/jperkel/gb_read/blob/main/circdna.gb), are included in this repository.)
5) By default, the program uses a three-letter amino acid code (e.g., Met for methionine). Use `-o` to use a one-character code instead, e.g., `cargo run -- -o -i pbr322.gb`.  
6) To create and view auto-generated documentation for this project, execute `cargo doc --open`.  
7) To run the program's test suite, execute `cargo test`.  

The R Markdown notebook, [Growth_of_cratesio.Rmd](https://github.com/jperkel/gb_read/blob/main/Growth_of_cratesio.Rmd), recreates the "[Rust rising](https://media.nature.com/lw800/magazine-assets/d41586-020-03382-2/d41586-020-03382-2_18629102.png)" graphic included in the article.
