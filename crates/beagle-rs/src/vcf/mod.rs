//! Port of the Java `vcf` package: the VCF data model (samples, markers, genotype
//! records), VCF/bref parsing, and VCF writing.
//!
//! This is the largest package; it is ported bottom-up across several chunks. Ported
//! so far: `Samples`, `Marker` (+ `MarkerUtils`, `MarkerParser`).

mod gt_rec;
mod marker;
mod marker_parser;
pub mod marker_utils;
mod samples;

pub use gt_rec::{allele_counts, allele_freq, to_vcf_rec, GTRec};
pub use marker::Marker;
pub use marker_parser::MarkerParser;
pub use samples::Samples;
