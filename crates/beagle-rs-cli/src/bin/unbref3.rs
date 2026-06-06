//! `unbref3` — expands a bref3 file back into VCF on standard output (drop-in replacement for
//! `java -jar unbref3.27Feb25.75f.jar`).
//!
//! Mirrors `bref/UnBref3.java`: arguments (excluding the program name) are forwarded verbatim to
//! the ported `UnBref3::main`.

use beagle_rs::bref::UnBref3;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    UnBref3::main(&args);
}
