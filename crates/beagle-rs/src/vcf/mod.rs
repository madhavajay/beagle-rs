//! Port of the Java `vcf` package: the VCF data model (samples, markers, genotype
//! records), VCF/bref parsing, and VCF writing.
//!
//! This is the largest package; it is ported bottom-up across several chunks. Ported
//! so far: `Samples`, `Marker` (+ `MarkerUtils`, `MarkerParser`).

mod marker;
mod marker_parser;
pub mod marker_utils;
mod samples;

pub use marker::Marker;
pub use marker_parser::MarkerParser;
pub use samples::Samples;
