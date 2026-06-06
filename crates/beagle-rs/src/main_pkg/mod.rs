//! Port of the Java `main` package — program constants, the `Par` parameter set, and the
//! top-level driver. (Named `main_pkg` because `main` is reserved for a binary entry point.)

mod par;
mod pedigree;
mod run_stats;
mod window_writer;

pub use par::{li_stephens_p_mismatch, Par};
pub use pedigree::Pedigree;
pub use run_stats::RunStats;
pub use window_writer::WindowWriter;

/// `Main.VERSION`.
pub const VERSION: &str = "(version 5.5)";
/// `Main.PROGRAM`.
pub const PROGRAM: &str = "beagle.27Feb25.75f.jar";
/// `Main.COMMAND`.
pub const COMMAND: &str = "java -jar beagle.27Feb25.75f.jar";
/// `Main.COPYRIGHT`.
pub const COPYRIGHT: &str = "Copyright (C) 2014-2024 Brian L. Browning";

/// `Main.SHORT_HELP` — the program name/version, copyright, and a one-line usage hint.
pub fn short_help() -> String {
    format!("{PROGRAM} {VERSION}\n{COPYRIGHT}\nEnter \"{COMMAND}\" to list command line argument")
}
