//! Command-line front-end for `beagle-rs`.
//!
//! Placeholder during the port: the CLI argument surface (`main/Par.java`) and the
//! phasing/imputation pipeline are ported in later phases. This binary currently
//! only reports build/version info so the workspace has a runnable target.

fn main() {
    eprintln!(
        "beagle-rs {} — Rust port of Beagle 5.5 (27Feb25.75f), in progress.",
        env!("CARGO_PKG_VERSION")
    );
    eprintln!("The CLI pipeline is not yet ported; see TODO.md for status.");
    std::process::exit(1);
}
