//! Port of `vcf/LowMafDiallelicGTRec.java` — the diallelic specialization of
//! `LowMafGTRec`, storing only the minor-allele haplotypes and missing samples.

use crate::ints::IntArray;

use super::{to_vcf_rec, GTRec, Marker, Samples};

/// Port of `vcf/LowMafDiallelicGTRec.java`.
pub struct LowMafDiallelicGTRec {
    marker: Marker,
    samples: Samples,
    n_haps: i32,
    major_allele: i32,
    minor_allele: i32,
    missing_samples: Vec<i32>,
    minor_alleles: Vec<i32>,
    is_phased: bool,
}

impl LowMafDiallelicGTRec {
    /// Builds a `LowMafDiallelicGTRec` from components (mirrors the `HapListRep` ctor).
    pub fn new(
        marker: Marker,
        samples: Samples,
        major_allele: i32,
        minor_alleles: Vec<i32>,
        missing_samples: Vec<i32>,
        is_phased: bool,
    ) -> Self {
        assert!(marker.n_alleles() == 2, "{}", marker.n_alleles());
        let n_haps = samples.size() << 1;
        LowMafDiallelicGTRec {
            marker,
            samples,
            n_haps,
            major_allele,
            minor_allele: 1 - major_allele,
            missing_samples,
            minor_alleles,
            is_phased,
        }
    }

    /// `majorAllele()`.
    pub fn major_allele(&self) -> i32 {
        self.major_allele
    }

    /// `alleleCount(int allele)`.
    pub fn allele_count(&self, allele: i32) -> i32 {
        if allele == self.major_allele {
            self.n_haps
                - self.minor_alleles.len() as i32
                - ((self.missing_samples.len() as i32) << 1)
        } else {
            self.minor_alleles.len() as i32
        }
    }
}

impl IntArray for LowMafDiallelicGTRec {
    fn size(&self) -> i32 {
        self.n_haps
    }

    fn get(&self, hap: i32) -> i32 {
        assert!(hap >= 0 && hap < self.n_haps, "{}", hap);
        if self.minor_alleles.binary_search(&hap).is_ok() {
            self.minor_allele
        } else if self.missing_samples.binary_search(&(hap >> 1)).is_ok() {
            -1
        } else {
            self.major_allele
        }
    }
}

impl GTRec for LowMafDiallelicGTRec {
    fn samples(&self) -> &Samples {
        &self.samples
    }
    fn marker(&self) -> &Marker {
        &self.marker
    }
    fn is_phased_sample(&self, sample: i32) -> bool {
        assert!(sample >= 0 && sample < self.samples.size(), "{}", sample);
        self.is_phased
    }
    fn is_phased(&self) -> bool {
        self.is_phased
    }
}

impl std::fmt::Display for LowMafDiallelicGTRec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&to_vcf_rec(self))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::MarkerParser;

    #[test]
    fn diallelic_low_maf_with_missing() {
        let marker = Marker::instance(
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/1",
            &MarkerParser::new(true, true, true, true),
        );
        let samples = Samples::new(
            &["s0".into(), "s1".into(), "s2".into()],
            &[true, true, true],
        );
        // major 0, minor 1 carried by haps [1,5]; sample 1 (haps 2,3) missing; unphased.
        let r = LowMafDiallelicGTRec::new(marker, samples, 0, vec![1, 5], vec![1], false);
        assert_eq!(r.size(), 6);
        assert_eq!(
            (0..6).map(|h| r.get(h)).collect::<Vec<_>>(),
            vec![0, 1, -1, -1, 0, 1]
        );
        assert_eq!(r.major_allele(), 0);
        assert_eq!(r.allele_count(1), 2);
        // major count = 6 - 2 minor - 2 missing = 2
        assert_eq!(r.allele_count(0), 2);
        assert!(!r.is_phased());
    }
}
