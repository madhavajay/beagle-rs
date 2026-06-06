//! `bref3` — converts a VCF file into bref3 format on standard output (drop-in replacement for
//! `java -jar bref3.27Feb25.75f.jar`).
//!
//! Mirrors `bref/Bref3.java`: arguments (excluding the program name) are forwarded verbatim to
//! the ported `Bref3::main`.

use beagle_rs::bref::Bref3;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    Bref3::main(&args);
}
