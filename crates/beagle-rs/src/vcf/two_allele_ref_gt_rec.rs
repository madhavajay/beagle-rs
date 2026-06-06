//! Port of `vcf/TwoAlleleRefGTRec.java` — an allele-coded reference record for a
//! diallelic marker, storing only the haplotypes that carry the minor allele.

use crate::ints::{packed_create, IndexArray, IntArray};

use super::allele_ref_gt_rec::check_indices_and_return_null_index;
use super::{to_vcf_rec, GTRec, Marker, RefGTRec, Samples};

/// Port of `vcf/TwoAlleleRefGTRec.java`.
pub struct TwoAlleleRefGTRec {
    marker: Marker,
    samples: Samples,
    n_haps: i32,
    major_allele: i32,
    minor_allele: i32,
    minor_alleles: Vec<i32>,
}

impl TwoAlleleRefGTRec {
    /// `new TwoAlleleRefGTRec(Marker, Samples, int[][] hapIndices)`.
    pub fn from_components(
        marker: Marker,
        samples: Samples,
        hap_indices: Vec<Option<Vec<i32>>>,
    ) -> Self {
        assert!(marker.n_alleles() == 2, "{}", marker.n_alleles());
        let n_haps = 2 * samples.size();
        let major_allele = check_indices_and_return_null_index(&hap_indices, n_haps);
        let minor_allele = 1 - major_allele;
        let minor_alleles = hap_indices[minor_allele as usize]
            .clone()
            .expect("minor allele row");
        TwoAlleleRefGTRec {
            marker,
            samples,
            n_haps,
            major_allele,
            minor_allele,
            minor_alleles,
        }
    }

    /// `new TwoAlleleRefGTRec(RefGTRec rec)`.
    pub fn from_ref_rec(rec: &dyn RefGTRec) -> Self {
        assert!(
            rec.marker().n_alleles() == 2,
            "{}",
            rec.marker().n_alleles()
        );
        let hap_indices = rec.allele_to_haps();
        let mut maj = 0usize;
        while hap_indices[maj].is_some() {
            maj += 1;
        }
        let major_allele = maj as i32;
        let minor_allele = 1 - major_allele;
        let minor_alleles = hap_indices[minor_allele as usize]
            .clone()
            .expect("minor allele row");
        TwoAlleleRefGTRec {
            marker: rec.marker().clone(),
            samples: rec.samples().clone(),
            n_haps: rec.size(),
            major_allele,
            minor_allele,
            minor_alleles,
        }
    }

    fn to_int_array(&self) -> Box<dyn IntArray> {
        let mut ia = vec![self.major_allele; self.n_haps as usize];
        for &h in &self.minor_alleles {
            ia[h as usize] = self.minor_allele;
        }
        packed_create(&ia, 2)
    }
}

impl IntArray for TwoAlleleRefGTRec {
    fn size(&self) -> i32 {
        self.n_haps
    }

    fn get(&self, hap: i32) -> i32 {
        assert!(hap >= 0 && hap < self.n_haps, "{}", hap);
        if self.minor_alleles.binary_search(&hap).is_ok() {
            self.minor_allele
        } else {
            self.major_allele
        }
    }
}

impl GTRec for TwoAlleleRefGTRec {
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

impl RefGTRec for TwoAlleleRefGTRec {
    fn allele_to_haps(&self) -> Vec<Option<Vec<i32>>> {
        let mut hap_indices: Vec<Option<Vec<i32>>> = vec![None, None];
        hap_indices[self.minor_allele as usize] = Some(self.minor_alleles.clone());
        hap_indices
    }

    fn hap_to_allele(&self) -> IndexArray {
        IndexArray::from_int_array(self.to_int_array(), 2)
    }

    fn n_allele_coded_haps(&self) -> i32 {
        self.minor_alleles.len() as i32
    }

    fn is_allele_coded(&self) -> bool {
        true
    }

    fn major_allele(&self) -> i32 {
        self.major_allele
    }

    fn allele_counts(&self) -> Vec<i32> {
        let mut cnts = vec![0i32; 2];
        cnts[self.major_allele as usize] = self.n_haps - self.minor_alleles.len() as i32;
        cnts[self.minor_allele as usize] = self.minor_alleles.len() as i32;
        cnts
    }

    fn allele_count(&self, allele: i32) -> i32 {
        if allele == self.major_allele {
            panic!("major allele");
        }
        self.minor_alleles.len() as i32
    }

    fn hap_index(&self, allele: i32, copy: i32) -> i32 {
        if allele == self.major_allele {
            panic!("major allele");
        }
        self.minor_alleles[copy as usize]
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

impl std::fmt::Display for TwoAlleleRefGTRec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&to_vcf_rec(self))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::MarkerParser;

    fn rec() -> TwoAlleleRefGTRec {
        let marker = Marker::instance(
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1",
            &MarkerParser::new(true, true, true, true),
        );
        let samples = Samples::new(&["s0".into(), "s1".into()], &[true, true]);
        TwoAlleleRefGTRec::from_components(marker, samples, vec![None, Some(vec![1, 3])])
    }

    #[test]
    fn diallelic_access() {
        let r = rec();
        assert_eq!(r.size(), 4);
        assert_eq!(r.major_allele(), 0);
        assert_eq!(
            (0..4).map(|h| r.get(h)).collect::<Vec<_>>(),
            vec![0, 1, 0, 1]
        );
        assert_eq!(r.allele_counts(), vec![2, 2]);
        assert_eq!(r.allele_count(1), 2);
        assert_eq!(r.hap_index(1, 1), 3);
        assert!(r.is_carrier(1, 1));
        assert_eq!(r.n_allele_coded_haps(), 2);
        let h2a = r.hap_to_allele();
        for h in 0..4 {
            assert_eq!(h2a.get(h), r.get(h));
        }
    }
}
