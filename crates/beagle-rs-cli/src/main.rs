//! `beagle-rs` — the Beagle phasing/imputation program (drop-in replacement for
//! `java -jar beagle.27Feb25.75f.jar`).
//!
//! Mirrors `main/Main.java`: command-line arguments (excluding the program name) are forwarded
//! verbatim to the ported driver, which parses them with `Par` and runs the pipeline.

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    beagle_rs::main_pkg::main(&args);
}
