//! `beagle-rs` — a Rust port of Beagle 5.5 (`27Feb25.75f`), genotype phasing and
//! genotype imputation.
//!
//! The port mirrors the original Java package structure one-to-one so each Rust
//! module maps directly to a Java package, and each type to a Java class. The goal
//! is byte-for-byte-identical output vs. the Java reference (see `TODO.md`).
//!
//! Ported packages (bottom-up):
//! - [`ints`] — integer-packed immutable arrays + maps (Java package `ints`).
//! - [`blbutil`] — base utilities: constants, float/double arrays, bit arrays, string
//!   splitting, BGZIP, file IO (Java package `blbutil`).

// The port mirrors Java's explicit bound checks (`x >= lo && x <= hi`) verbatim so the
// Rust lines correspond directly to the original source; clippy reads that as a
// `manual_range_contains` candidate. Keeping the Java form aids review/parity.
#![allow(clippy::manual_range_contains)]

pub mod beagleutil;
pub mod blbutil;
pub mod ints;
pub mod jdk;
