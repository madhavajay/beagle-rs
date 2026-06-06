//! Port of the Java `vcf` package: the VCF data model (samples, markers, genotype
//! records), VCF/bref parsing, and VCF writing.
//!
//! This is the largest package; it is ported bottom-up across several chunks. Ported
//! so far: `Samples`.

mod samples;

pub use samples::Samples;
