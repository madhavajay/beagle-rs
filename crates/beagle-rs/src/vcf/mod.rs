//! Port of the Java `vcf` package: the VCF data model (samples, markers, genotype
//! records), VCF/bref parsing, and VCF writing.
//!
//! This is the largest package; it is ported bottom-up across several chunks. Ported
//! so far: `Samples`, `Marker` (+ `MarkerUtils`, `MarkerParser`).

mod allele_ref_gt_rec;
mod basic_gt_rec;
mod bit_array_ref_gt_rec;
mod gt;
mod gt_rec;
mod hap_ref_gt_rec;
mod int_array_ref_gt_rec;
mod low_maf_diallelic_gt_rec;
mod low_maf_gt_rec;
mod marker;
mod marker_parser;
pub mod marker_utils;
mod markers;
mod ref_gt_rec;
mod samples;
mod two_allele_ref_gt_rec;
mod vcf_header;
mod vcf_meta_info;

pub use allele_ref_gt_rec::AlleleRefGTRec;
pub use basic_gt_rec::BasicGTRec;
pub use bit_array_ref_gt_rec::BitArrayRefGTRec;
pub use gt::GT;
pub use gt_rec::{allele_counts, allele_freq, to_vcf_rec, GTRec};
pub use hap_ref_gt_rec::HapRefGTRec;
pub use int_array_ref_gt_rec::IntArrayRefGTRec;
pub use low_maf_diallelic_gt_rec::LowMafDiallelicGTRec;
pub use low_maf_gt_rec::LowMafGTRec;
pub use marker::Marker;
pub use marker_parser::MarkerParser;
pub use markers::Markers;
pub use ref_gt_rec::{allele_ref_gt_rec_from_components, allele_ref_gt_rec_from_rec, RefGTRec};
pub use samples::Samples;
pub use two_allele_ref_gt_rec::TwoAlleleRefGTRec;
pub use vcf_header::{VcfHeader, HEADER_PREFIX};
pub use vcf_meta_info::VcfMetaInfo;
