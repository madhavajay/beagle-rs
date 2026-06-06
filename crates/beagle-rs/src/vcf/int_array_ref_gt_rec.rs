//! Port of `vcf/IntArrayRefGTRec.java` — a map-coded (NOT allele-coded) reference
//! record storing the allele on each haplotype.
//!
//! Java stores the alleles as a single packed `IntArray`. Here we store the unpacked
//! alleles and reconstruct the packed array on demand for `maps()`/`map()`/
//! `hapToAllele()`. `get()` and all returned values are identical to Java (and the
//! reconstructed packed array is byte-identical to Java's stored one).

use crate::ints::{packed_create, IndexArray, IntArray};

use super::ref_gt_rec::non_null_cnt;
use super::{to_vcf_rec, GTRec, Marker, RefGTRec, Samples};

/// Port of `vcf/IntArrayRefGTRec.java`.
pub struct IntArrayRefGTRec {
    marker: Marker,
    samples: Samples,
    n_alleles: i32,
    alleles: Vec<i32>,
}

fn major_allele_of(al_cnts: &[i32]) -> i32 {
    let mut maj = 0usize;
    for al in 1..al_cnts.len() {
        if al_cnts[al] > al_cnts[maj] {
            maj = al;
        }
    }
    maj as i32
}

impl IntArrayRefGTRec {
    /// `new IntArrayRefGTRec(Marker, Samples, int[] alleles)`.
    pub fn from_alleles(marker: Marker, samples: Samples, alleles: &[i32]) -> Self {
        assert!(
            alleles.len() == 2 * samples.size() as usize,
            "{}",
            alleles.len()
        );
        let n_alleles = marker.n_alleles();
        for &a in alleles {
            assert!(a >= 0 && a < n_alleles, "{}", a);
        }
        IntArrayRefGTRec {
            marker,
            samples,
            n_alleles,
            alleles: alleles.to_vec(),
        }
    }

    /// `new IntArrayRefGTRec(Marker, Samples, IndexArray alleles)`.
    pub fn from_index_array(marker: Marker, samples: Samples, alleles: &IndexArray) -> Self {
        assert!(alleles.size() == 2 * samples.size(), "{}", alleles.size());
        assert!(
            alleles.value_size() <= marker.n_alleles(),
            "{}",
            alleles.value_size()
        );
        let n_alleles = marker.n_alleles();
        let raw: Vec<i32> = (0..alleles.size()).map(|h| alleles.get(h)).collect();
        IntArrayRefGTRec {
            marker,
            samples,
            n_alleles,
            alleles: raw,
        }
    }

    fn packed(&self) -> Box<dyn IntArray> {
        packed_create(&self.alleles, self.n_alleles)
    }
}

impl IntArray for IntArrayRefGTRec {
    fn size(&self) -> i32 {
        2 * self.samples.size()
    }

    fn get(&self, hap: i32) -> i32 {
        self.alleles[hap as usize]
    }
}

impl GTRec for IntArrayRefGTRec {
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

impl RefGTRec for IntArrayRefGTRec {
    fn allele_to_haps(&self) -> Vec<Option<Vec<i32>>> {
        let al_cnts = self.allele_counts();
        let maj_allele = major_allele_of(&al_cnts) as usize;
        let mut hap_indices: Vec<Option<Vec<i32>>> = (0..al_cnts.len())
            .map(|al| {
                if al != maj_allele {
                    Some(Vec::with_capacity(al_cnts[al] as usize))
                } else {
                    None
                }
            })
            .collect();
        for (h, &al) in self.alleles.iter().enumerate() {
            let al = al as usize;
            if al != maj_allele {
                hap_indices[al].as_mut().unwrap().push(h as i32);
            }
        }
        hap_indices
    }

    fn hap_to_allele(&self) -> IndexArray {
        IndexArray::from_int_array(self.packed(), self.n_alleles)
    }

    fn n_allele_coded_haps(&self) -> i32 {
        non_null_cnt(&self.allele_to_haps())
    }

    fn is_allele_coded(&self) -> bool {
        false
    }

    fn major_allele(&self) -> i32 {
        major_allele_of(&self.allele_counts())
    }

    fn allele_counts(&self) -> Vec<i32> {
        let mut al_cnts = vec![0i32; self.n_alleles as usize];
        for &al in &self.alleles {
            al_cnts[al as usize] += 1;
        }
        al_cnts
    }

    fn allele_count(&self, allele: i32) -> i32 {
        let al_cnts = self.allele_counts();
        if allele == major_allele_of(&al_cnts) {
            panic!("major allele");
        }
        al_cnts[allele as usize]
    }

    fn hap_index(&self, allele: i32, copy: i32) -> i32 {
        let hap_indices = self.allele_to_haps();
        match &hap_indices[allele as usize] {
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
        vec![self.packed()]
    }

    fn map(&self, index: i32) -> Box<dyn IntArray> {
        assert!(index == 0, "{}", index);
        self.packed()
    }
}

impl std::fmt::Display for IntArrayRefGTRec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&to_vcf_rec(self))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::MarkerParser;

    fn rec() -> IntArrayRefGTRec {
        // triallelic; alleles per hap: [0,1,0,2,0,0] -> major 0
        let marker = Marker::instance(
            "chr1\t100\t.\tA\tC,G\t.\tPASS\t.\tGT\t0|1",
            &MarkerParser::new(true, true, true, true),
        );
        let samples = Samples::new(
            &["s0".into(), "s1".into(), "s2".into()],
            &[true, true, true],
        );
        IntArrayRefGTRec::from_alleles(marker, samples, &[0, 1, 0, 2, 0, 0])
    }

    #[test]
    fn map_coded_access() {
        let r = rec();
        assert_eq!(r.size(), 6);
        assert!(!r.is_allele_coded());
        assert_eq!(r.get(1), 1);
        assert_eq!(r.get(3), 2);
        assert_eq!(r.allele_counts(), vec![4, 1, 1]);
        assert_eq!(r.major_allele(), 0);
        assert_eq!(r.allele_count(1), 1);
        assert_eq!(r.hap_index(2, 0), 3);
        assert!(r.is_carrier(2, 3));
        assert_eq!(r.n_allele_coded_haps(), 2);
        let h2a = r.hap_to_allele();
        for h in 0..6 {
            assert_eq!(h2a.get(h), r.get(h));
        }
        let m = r.map(0);
        for h in 0..6 {
            assert_eq!(m.get(h), r.get(h));
        }
        // allele_to_haps: major row None, allele1=[1], allele2=[3]
        let a2h = r.allele_to_haps();
        assert_eq!(a2h[0], None);
        assert_eq!(a2h[1], Some(vec![1]));
        assert_eq!(a2h[2], Some(vec![3]));
    }
}
