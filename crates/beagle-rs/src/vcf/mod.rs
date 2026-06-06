//! Port of the Java `vcf` package: the VCF data model (samples, markers, genotype
//! records), VCF/bref parsing, and VCF writing.
//!
//! This is the largest package; it is ported bottom-up across several chunks. Ported
//! so far: `Samples`, `Marker` (+ `MarkerUtils`, `MarkerParser`).

mod allele_ref_gt_rec;
mod basic_gt;
mod basic_gt_rec;
mod bit_array_gt_rec;
mod bit_array_ref_gt_rec;
mod genetic_map;
mod gt;
mod gt_rec;
mod hap_ref_gt_rec;
mod int_array_ref_gt_rec;
mod low_maf_diallelic_gt_rec;
mod low_maf_gt_rec;
mod marker;
mod marker_indices;
mod marker_map;
mod marker_parser;
pub mod marker_utils;
mod markers;
mod plink_gen_map;
mod position_map;
mod ref_gt;
mod ref_gt_rec;
mod restricted_gt;
mod samples;
mod spliced_gt;
mod steps;
mod two_allele_ref_gt_rec;
mod vcf_header;
mod vcf_meta_info;
mod vcf_rec;
mod vcf_rec_builder;
mod vcf_rec_gt_parser;
mod vcf_writer;
mod xref_gt;

pub use allele_ref_gt_rec::AlleleRefGTRec;
pub use basic_gt::{genotype, BasicGT};
pub use basic_gt_rec::BasicGTRec;
pub use bit_array_gt_rec::BitArrayGTRec;
pub use bit_array_ref_gt_rec::BitArrayRefGTRec;
pub use genetic_map::{
    gen_pos_markers, gen_pos_markers_min_dist, genetic_map_from_file, GeneticMap,
};
pub use gt::GT;
pub use gt_rec::{allele_counts, allele_freq, to_vcf_rec, GTRec};
pub use hap_ref_gt_rec::HapRefGTRec;
pub use int_array_ref_gt_rec::IntArrayRefGTRec;
pub use low_maf_diallelic_gt_rec::LowMafDiallelicGTRec;
pub use low_maf_gt_rec::LowMafGTRec;
pub use marker::Marker;
pub use marker_indices::MarkerIndices;
pub use marker_map::{mean_single_base_gen_dist, MarkerMap};
pub use marker_parser::MarkerParser;
pub use markers::Markers;
pub use plink_gen_map::PlinkGenMap;
pub use position_map::PositionMap;
pub use ref_gt::RefGT;
pub use ref_gt_rec::{allele_ref_gt_rec_from_components, allele_ref_gt_rec_from_rec, RefGTRec};
pub use restricted_gt::RestrictedGT;
pub use samples::Samples;
pub use spliced_gt::SplicedGT;
pub use steps::Steps;
pub use two_allele_ref_gt_rec::TwoAlleleRefGTRec;
pub use vcf_header::{VcfHeader, HEADER_PREFIX};
pub use vcf_meta_info::VcfMetaInfo;
pub use vcf_rec::{gt_index, VcfRec};
pub use vcf_rec_builder::VcfRecBuilder;
pub use vcf_rec_gt_parser::{HapListRep, VcfRecGTParser};
pub use vcf_writer::{
    append_records_gt, append_records_recs, print_fixed_fields_gt, write_meta_lines,
    write_meta_lines_gt,
};
pub use xref_gt::XRefGT;
