//! Port of `vcf/RefGTRec.java` — phased, non-missing genotypes for one marker, with two
//! representation modes (allele-coded: store haplotypes per non-major allele; or
//! map-coded). Implementations: `AlleleRefGTRec`, `TwoAlleleRefGTRec`, `BitArrayRefGTRec`,
//! `IntArrayRefGTRec`, `HapRefGTRec`, `LowMafGTRec`, `LowMafDiallelicGTRec` (ported across
//! chunks). The static `alleleRefGTRec` factories are added once the impls exist.

use crate::ints::{IndexArray, IntArray};

use super::GTRec;

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
}

/// Sum of the lengths of the non-`None` rows (number of non-major-allele haplotypes).
pub(crate) fn non_null_cnt(allele_to_haps: &[Option<Vec<i32>>]) -> i32 {
    allele_to_haps
        .iter()
        .filter_map(|r| r.as_ref().map(|v| v.len() as i32))
        .sum()
}
