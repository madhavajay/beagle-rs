//! Port of `vcf/LowMafGTRec.java` — genotypes for one marker storing, per non-major
//! allele, the haplotypes carrying it, plus the samples with missing genotypes. All
//! genotypes are treated as unphased if any sample is unphased/missing.
//!
//! The `VcfRecGTParser.HapListRep` constructor is added with the parser; the component
//! constructor is provided here.

use crate::ints::IntArray;

use super::{to_vcf_rec, GTRec, Marker, Samples};

/// Port of `vcf/LowMafGTRec.java`.
pub struct LowMafGTRec {
    marker: Marker,
    samples: Samples,
    n_haps: i32,
    major_allele: i32,
    hap_indices: Vec<Option<Vec<i32>>>, // major-allele row is None
    missing_samples: Vec<i32>,
    is_phased: bool,
}

impl LowMafGTRec {
    /// Builds a `LowMafGTRec` from components (mirrors the `HapListRep` constructor).
    pub fn new(
        marker: Marker,
        samples: Samples,
        major_allele: i32,
        hap_indices: Vec<Option<Vec<i32>>>,
        missing_samples: Vec<i32>,
        is_phased: bool,
    ) -> Self {
        let n_haps = samples.size() << 1;
        LowMafGTRec {
            marker,
            samples,
            n_haps,
            major_allele,
            hap_indices,
            missing_samples,
            is_phased,
        }
    }

    /// `new LowMafGTRec(VcfRecGTParser.HapListRep listRep)`.
    pub fn from_hap_list_rep(list_rep: &super::HapListRep) -> Self {
        LowMafGTRec::new(
            list_rep.marker().clone(),
            list_rep.samples().clone(),
            list_rep.major_allele(),
            list_rep.hap_lists(true),
            list_rep.missing_samples(),
            list_rep.is_phased(),
        )
    }

    /// `majorAllele()`.
    pub fn major_allele(&self) -> i32 {
        self.major_allele
    }

    /// `alleleCount(int allele)`. NOTE (vcf-2): the major-allele branch preserves a Java
    /// bug — it subtracts `hapIndices.length` (the allele count) per non-major allele
    /// instead of that row's length, so the major-allele count is generally wrong.
    pub fn allele_count(&self, allele: i32) -> i32 {
        if allele == self.major_allele {
            let mut n = self.n_haps - ((self.missing_samples.len() as i32) << 1);
            for al in 0..self.hap_indices.len() {
                if al as i32 != self.major_allele {
                    n -= self.hap_indices.len() as i32; // vcf-2: Java bug, preserved
                }
            }
            n
        } else {
            self.hap_indices[allele as usize]
                .as_ref()
                .expect("non-major allele row")
                .len() as i32
        }
    }
}

impl IntArray for LowMafGTRec {
    fn size(&self) -> i32 {
        self.n_haps
    }

    fn get(&self, hap: i32) -> i32 {
        assert!(hap >= 0 && hap < self.n_haps, "{}", hap);
        for (j, row) in self.hap_indices.iter().enumerate() {
            if j as i32 != self.major_allele {
                if let Some(ia) = row {
                    if ia.binary_search(&hap).is_ok() {
                        return j as i32;
                    }
                }
            }
        }
        if self.missing_samples.binary_search(&(hap >> 1)).is_ok() {
            -1
        } else {
            self.major_allele
        }
    }
}

impl GTRec for LowMafGTRec {
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

impl std::fmt::Display for LowMafGTRec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&to_vcf_rec(self))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::MarkerParser;

    #[test]
    fn low_maf_get_with_missing() {
        // triallelic; major allele 0. allele1 haps [1], allele2 haps [5]; sample 1
        // (haps 2,3) missing; unphased.
        let marker = Marker::instance(
            "chr1\t100\t.\tA\tC,G\t.\tPASS\t.\tGT\t0/1",
            &MarkerParser::new(true, true, true, true),
        );
        let samples = Samples::new(
            &["s0".into(), "s1".into(), "s2".into()],
            &[true, true, true],
        );
        let r = LowMafGTRec::new(
            marker,
            samples,
            0,
            vec![None, Some(vec![1]), Some(vec![5])],
            vec![1], // sample 1 missing
            false,
        );
        assert_eq!(r.size(), 6);
        assert_eq!(r.get(0), 0);
        assert_eq!(r.get(1), 1);
        assert_eq!(r.get(2), -1); // sample 1 missing
        assert_eq!(r.get(3), -1);
        assert_eq!(r.get(4), 0);
        assert_eq!(r.get(5), 2);
        assert!(!r.is_phased());
        assert_eq!(r.major_allele(), 0);
        assert_eq!(r.allele_count(1), 1);
    }
}
