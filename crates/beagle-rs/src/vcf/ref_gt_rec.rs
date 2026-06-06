//! Port of `vcf/RefGTRec.java` — phased, non-missing genotypes for one marker, with two
//! representation modes (allele-coded: store haplotypes per non-major allele; or
//! map-coded). Implementations: `AlleleRefGTRec`, `TwoAlleleRefGTRec`, `BitArrayRefGTRec`,
//! `IntArrayRefGTRec`, `HapRefGTRec`, `LowMafGTRec`, `LowMafDiallelicGTRec` (ported across
//! chunks). The static `alleleRefGTRec` factories are added once the impls exist.

use crate::ints::{IndexArray, IntArray};

use super::{AlleleRefGTRec, GTRec, Marker, Samples, TwoAlleleRefGTRec};

/// Port of the `vcf/RefGTRec.java` interface.
pub trait RefGTRec: GTRec {
    /// `alleleToHaps()` — element `j` is the sorted haplotypes carrying allele `j`, or
    /// `None` for the major allele (exactly one `None`).
    fn allele_to_haps(&self) -> Vec<Option<Vec<i32>>>;

    /// `hapToAllele()` — maps haplotype → allele.
    fn hap_to_allele(&self) -> IndexArray;

    /// `nAlleleCodedHaps()` — sum of lengths of non-`None` rows of `alleleToHaps()`.
    fn n_allele_coded_haps(&self) -> i32;

    /// `isAlleleCoded()`.
    fn is_allele_coded(&self) -> bool;

    /// `majorAllele()`.
    fn major_allele(&self) -> i32;

    /// `alleleCounts()`.
    fn allele_counts(&self) -> Vec<i32>;

    /// `alleleCount(int allele)` — count for a non-major allele.
    fn allele_count(&self, allele: i32) -> i32;

    /// `hapIndex(int allele, int copy)`.
    fn hap_index(&self, allele: i32, copy: i32) -> i32;

    /// `isCarrier(int allele, int hap)`.
    fn is_carrier(&self, allele: i32, hap: i32) -> bool;

    /// `nMaps()`.
    fn n_maps(&self) -> i32;

    /// `maps()` — composed maps from haplotype index to allele.
    fn maps(&self) -> Vec<Box<dyn IntArray>>;

    /// `map(int index)`.
    fn map(&self, index: i32) -> Box<dyn IntArray>;

    /// Identity key for a sequence-coded record's `hapToSeq` map (the `Rc`'s address), used by
    /// the bref3 writer to group records sharing one map into a block. `None` for allele-coded
    /// records (Java compares the `IntArray` object reference from `map(0)`).
    fn seq_block_key(&self) -> Option<usize> {
        None
    }
}

/// `RefGTRec.alleleRefGTRec(Marker, Samples, int[][])` — biallelic → `TwoAlleleRefGTRec`,
/// otherwise `AlleleRefGTRec`.
pub fn allele_ref_gt_rec_from_components(
    marker: Marker,
    samples: Samples,
    allele_to_haps: Vec<Option<Vec<i32>>>,
) -> Box<dyn RefGTRec> {
    if marker.n_alleles() == 2 {
        Box::new(TwoAlleleRefGTRec::from_components(
            marker,
            samples,
            allele_to_haps,
        ))
    } else {
        Box::new(AlleleRefGTRec::from_components(
            marker,
            samples,
            allele_to_haps,
        ))
    }
}

/// `RefGTRec.alleleRefGTRec(VcfRecGTParser)` — parses a VCF record's GT field into an
/// allele-coded record (biallelic → `TwoAlleleRefGTRec`, else `AlleleRefGTRec`).
pub fn allele_ref_gt_rec_from_parser(gtp: &super::VcfRecGTParser) -> Box<dyn RefGTRec> {
    let marker = gtp.marker().clone();
    let samples = gtp.samples().clone();
    let non_maj = gtp.non_maj_ref_indices();
    if gtp.n_alleles() == 2 {
        Box::new(TwoAlleleRefGTRec::from_components(marker, samples, non_maj))
    } else {
        Box::new(AlleleRefGTRec::from_components(marker, samples, non_maj))
    }
}

/// `RefGTRec.alleleRefGTRec(RefGTRec)` — returns an allele-coded record. (Unlike Java,
/// this always reconstructs rather than returning an already-allele-coded `rec` as-is;
/// the result is identical for the immutable record.)
pub fn allele_ref_gt_rec_from_rec(rec: &dyn RefGTRec) -> Box<dyn RefGTRec> {
    if rec.marker().n_alleles() == 2 {
        Box::new(TwoAlleleRefGTRec::from_ref_rec(rec))
    } else {
        Box::new(AlleleRefGTRec::from_ref_rec(rec))
    }
}

/// Sum of the lengths of the non-`None` rows (number of non-major-allele haplotypes).
pub(crate) fn non_null_cnt(allele_to_haps: &[Option<Vec<i32>>]) -> i32 {
    allele_to_haps
        .iter()
        .filter_map(|r| r.as_ref().map(|v| v.len() as i32))
        .sum()
}
