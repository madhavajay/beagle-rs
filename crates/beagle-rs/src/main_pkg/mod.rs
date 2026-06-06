//! Port of the Java `main` package — program constants, the `Par` parameter set, and the
//! top-level driver. (Named `main_pkg` because `main` is reserved for a binary entry point.)

mod par;

pub use par::{li_stephens_p_mismatch, Par};

/// `Main.VERSION`.
pub const VERSION: &str = "(version 5.5)";
/// `Main.PROGRAM`.
pub const PROGRAM: &str = "beagle.27Feb25.75f.jar";
/// `Main.COMMAND`.
pub const COMMAND: &str = "java -jar beagle.27Feb25.75f.jar";
/// `Main.COPYRIGHT`.
pub const COPYRIGHT: &str = "Copyright (C) 2014-2024 Brian L. Browning";
