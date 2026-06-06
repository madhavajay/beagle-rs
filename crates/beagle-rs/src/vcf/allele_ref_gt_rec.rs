//! Port of `vcf/AlleleRefGTRec.java` — an allele-coded reference record storing, for
//! each non-major allele, the sorted haplotypes carrying it.
//!
//! The `VcfRecGTParser`-based constructor is added with that parser. The
//! component constructor and the `RefGTRec`-copy constructor are ported here.

use crate::ints::{packed_create, IndexArray, IntArray};

use super::ref_gt_rec::non_null_cnt;
use super::{to_vcf_rec, GTRec, Marker, RefGTRec, Samples};

/// Port of `vcf/AlleleRefGTRec.java`.
pub struct AlleleRefGTRec {
    marker: Marker,
    samples: Samples,
    n_haps: i32,
    major_allele: i32,
    allele_to_haps: Vec<Option<Vec<i32>>>,
}

fn check_sorted(ia: &[i32], n_haps: i32) {
    if !ia.is_empty() && (ia[0] < 0 || ia[ia.len() - 1] >= n_haps) {
        panic!("invalid array");
    }
    for k in 1..ia.len() {
        if ia[k - 1] >= ia[k] {
            panic!("invalid array");
        }
    }
}

/// `checkIndicesAndReturnNullIndex` — validates rows and returns the single `None` index.
fn check_indices_and_return_null_index(hap_indices: &[Option<Vec<i32>>], n_haps: i32) -> i32 {
    let mut maj_allele: i32 = -1;
    for (j, row) in hap_indices.iter().enumerate() {
        match row {
            None => {
                if maj_allele == -1 {
                    maj_allele = j as i32;
                } else {
                    panic!("invalid array");
                }
            }
            Some(ia) => check_sorted(ia, n_haps),
        }
    }
    if maj_allele == -1 {
        panic!("invalid array");
    }
    maj_allele
}

impl AlleleRefGTRec {
    /// `new AlleleRefGTRec(Marker, Samples, int[][] hapIndices)`.
    pub fn from_components(
        marker: Marker,
        samples: Samples,
        hap_indices: Vec<Option<Vec<i32>>>,
    ) -> Self {
        let n_haps = 2 * samples.size();
        let major_allele = check_indices_and_return_null_index(&hap_indices, n_haps);
        AlleleRefGTRec {
            marker,
            samples,
            n_haps,
            major_allele,
            allele_to_haps: hap_indices,
        }
    }

    /// `new AlleleRefGTRec(RefGTRec rec)`.
    pub fn from_ref_rec(rec: &dyn RefGTRec) -> Self {
        let allele_to_haps = rec.allele_to_haps();
        let mut maj = 0usize;
        while allele_to_haps[maj].is_some() {
            maj += 1;
        }
        AlleleRefGTRec {
            marker: rec.marker().clone(),
            samples: rec.samples().clone(),
            n_haps: rec.size(),
            major_allele: maj as i32,
            allele_to_haps,
        }
    }

    fn to_int_array(&self) -> Box<dyn IntArray> {
        let mut ia = vec![self.major_allele; self.n_haps as usize];
        for (al, row) in self.allele_to_haps.iter().enumerate() {
            if let Some(haps) = row {
                for &h in haps {
                    ia[h as usize] = al as i32;
                }
            }
        }
        packed_create(&ia, self.allele_to_haps.len() as i32)
    }
}

impl IntArray for AlleleRefGTRec {
    fn size(&self) -> i32 {
        self.n_haps
    }

    fn get(&self, hap: i32) -> i32 {
        assert!(hap >= 0 && hap < self.n_haps, "{}", hap);
        for (j, row) in self.allele_to_haps.iter().enumerate() {
            if j as i32 != self.major_allele {
                if let Some(ia) = row {
                    if ia.binary_search(&hap).is_ok() {
                        return j as i32;
                    }
                }
            }
        }
        self.major_allele
    }
}

impl GTRec for AlleleRefGTRec {
    fn samples(&self) -> &Samples {
        &self.samples
    }
    fn marker(&self) -> &Marker {
        &self.marker
    }
    fn is_phased_sample(&self, sample: i32) -> bool {
        assert!(sample >= 0 && sample < self.samples.size(), "{}", sample);
        true
    }
    fn is_phased(&self) -> bool {
        true
    }
}

impl RefGTRec for AlleleRefGTRec {
    fn allele_to_haps(&self) -> Vec<Option<Vec<i32>>> {
        self.allele_to_haps.clone()
    }

    fn hap_to_allele(&self) -> IndexArray {
        IndexArray::from_int_array(self.to_int_array(), self.allele_to_haps.len() as i32)
    }

    fn n_allele_coded_haps(&self) -> i32 {
        non_null_cnt(&self.allele_to_haps)
    }

    fn is_allele_coded(&self) -> bool {
        true
    }

    fn major_allele(&self) -> i32 {
        self.major_allele
    }

    fn allele_counts(&self) -> Vec<i32> {
        let mut cnts = vec![0i32; self.marker.n_alleles() as usize];
        cnts[self.major_allele as usize] = self.n_haps;
        for (j, row) in self.allele_to_haps.iter().enumerate() {
            if j as i32 != self.major_allele {
                let cnt = row.as_ref().expect("non-major allele row").len() as i32;
                cnts[j] = cnt;
                cnts[self.major_allele as usize] -= cnt;
            }
        }
        cnts
    }

    fn allele_count(&self, allele: i32) -> i32 {
        match &self.allele_to_haps[allele as usize] {
            None => panic!("major allele"),
            Some(ia) => ia.len() as i32,
        }
    }

    fn hap_index(&self, allele: i32, copy: i32) -> i32 {
        match &self.allele_to_haps[allele as usize] {
            None => panic!("major allele"),
            Some(ia) => ia[copy as usize],
        }
    }

    fn is_carrier(&self, allele: i32, hap: i32) -> bool {
        self.get(hap) == allele
    }

    fn n_maps(&self) -> i32 {
        1
    }

    fn maps(&self) -> Vec<Box<dyn IntArray>> {
        vec![self.to_int_array()]
    }

    fn map(&self, index: i32) -> Box<dyn IntArray> {
        assert!(index == 0, "{}", index);
        self.to_int_array()
    }
}

impl std::fmt::Display for AlleleRefGTRec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&to_vcf_rec(self))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::{allele_counts, MarkerParser};

    fn rec() -> AlleleRefGTRec {
        // 4 haplotypes; biallelic marker. major allele 0; allele 1 carried by haps 1,3.
        let marker = Marker::instance(
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1",
            &MarkerParser::new(true, true, true, true),
        );
        let samples = Samples::new(&["s0".into(), "s1".into()], &[true, true]);
        AlleleRefGTRec::from_components(marker, samples, vec![None, Some(vec![1, 3])])
    }

    #[test]
    fn allele_coded_access() {
        let r = rec();
        assert_eq!(r.size(), 4);
        assert_eq!(r.major_allele(), 0);
        assert_eq!(r.get(0), 0);
        assert_eq!(r.get(1), 1);
        assert_eq!(r.get(2), 0);
        assert_eq!(r.get(3), 1);
        assert!(r.is_allele_coded());
        assert_eq!(r.allele_counts(), vec![2, 2]);
        assert_eq!(r.allele_count(1), 2);
        assert_eq!(r.hap_index(1, 0), 1);
        assert_eq!(r.hap_index(1, 1), 3);
        assert!(r.is_carrier(1, 3));
        assert!(!r.is_carrier(1, 0));
        assert_eq!(r.n_allele_coded_haps(), 2);
        // hapToAllele round-trips via get
        let h2a = r.hap_to_allele();
        for h in 0..4 {
            assert_eq!(h2a.get(h), r.get(h));
        }
        // alleleCounts via the generic GTRec helper agrees
        assert_eq!(allele_counts(&r), vec![2, 2]);
    }

    #[test]
    fn maps_compose_to_alleles() {
        let r = rec();
        assert_eq!(r.n_maps(), 1);
        let m = r.map(0);
        for h in 0..4 {
            assert_eq!(m.get(h), r.get(h));
        }
    }
}
