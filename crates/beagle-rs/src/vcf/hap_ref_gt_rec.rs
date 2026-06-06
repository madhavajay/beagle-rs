//! Port of `vcf/HapRefGTRec.java` — a sequence-coded reference record: each haplotype
//! maps to a distinct allele-sequence index (`hapToSeq`), and each sequence maps to the
//! marker allele (`seqToAllele`). `get(h) = seqToAllele[hapToSeq[h]]`.
//!
//! As with `IntArrayRefGTRec`, the two maps are stored unpacked and the packed `IntArray`s
//! are reconstructed on demand for `maps()`/`map()` (value-identical to Java).

use std::rc::Rc;

use crate::ints::{packed_create, IndexArray, IntArray};

use super::ref_gt_rec::non_null_cnt;
use super::{to_vcf_rec, GTRec, Marker, RefGTRec, Samples};

/// Port of `vcf/HapRefGTRec.java`.
///
/// `hap_to_seq` is held in an `Rc` so that all records produced by a single
/// `SeqCoder3::get_compressed_list` group share one map by identity — the bref3 writer groups
/// records into a block by comparing this identity (Java compares the `IntArray` object).
pub struct HapRefGTRec {
    marker: Marker,
    samples: Samples,
    hap_to_seq: Rc<Vec<i32>>,
    seq_to_allele: Vec<i32>,
    n_alleles: i32,
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

impl HapRefGTRec {
    /// `new HapRefGTRec(Marker, Samples, IntArray hapToSeq, IntArray seqToAllele)`.
    pub fn new(
        marker: Marker,
        samples: Samples,
        hap_to_seq: &dyn IntArray,
        seq_to_allele: &dyn IntArray,
    ) -> Self {
        let hts: Vec<i32> = (0..hap_to_seq.size()).map(|i| hap_to_seq.get(i)).collect();
        Self::new_shared(marker, samples, Rc::new(hts), seq_to_allele)
    }

    /// Like [`HapRefGTRec::new`] but shares an existing `hap_to_seq` map (`Rc`) so that records
    /// in the same compressed group compare equal by identity for bref3 block grouping.
    pub fn new_shared(
        marker: Marker,
        samples: Samples,
        hap_to_seq: Rc<Vec<i32>>,
        seq_to_allele: &dyn IntArray,
    ) -> Self {
        assert!(
            hap_to_seq.len() as i32 == 2 * samples.size(),
            "inconsistent data"
        );
        let n_alleles = marker.n_alleles();
        let sta: Vec<i32> = (0..seq_to_allele.size())
            .map(|i| seq_to_allele.get(i))
            .collect();
        HapRefGTRec {
            marker,
            samples,
            hap_to_seq,
            seq_to_allele: sta,
            n_alleles,
        }
    }

    fn n_seqs(&self) -> i32 {
        self.seq_to_allele.len() as i32
    }

    fn map0(&self) -> Box<dyn IntArray> {
        packed_create(&self.hap_to_seq, self.n_seqs())
    }

    fn map1(&self) -> Box<dyn IntArray> {
        packed_create(&self.seq_to_allele, self.n_alleles)
    }

    fn alleles_vec(&self) -> Vec<i32> {
        self.hap_to_seq
            .iter()
            .map(|&seq| self.seq_to_allele[seq as usize])
            .collect()
    }
}

impl IntArray for HapRefGTRec {
    fn size(&self) -> i32 {
        self.hap_to_seq.len() as i32
    }

    fn get(&self, hap: i32) -> i32 {
        self.seq_to_allele[self.hap_to_seq[hap as usize] as usize]
    }
}

impl GTRec for HapRefGTRec {
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

impl RefGTRec for HapRefGTRec {
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
        for h in 0..self.size() {
            let al = self.get(h) as usize;
            if al != maj_allele {
                hap_indices[al].as_mut().unwrap().push(h);
            }
        }
        hap_indices
    }

    fn hap_to_allele(&self) -> IndexArray {
        IndexArray::from_int_array(
            packed_create(&self.alleles_vec(), self.n_alleles),
            self.n_alleles,
        )
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
        for h in 0..self.size() {
            al_cnts[self.get(h) as usize] += 1;
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
        2
    }

    fn maps(&self) -> Vec<Box<dyn IntArray>> {
        vec![self.map0(), self.map1()]
    }

    fn map(&self, index: i32) -> Box<dyn IntArray> {
        match index {
            0 => self.map0(),
            1 => self.map1(),
            _ => panic!("{}", index),
        }
    }

    fn seq_block_key(&self) -> Option<usize> {
        Some(Rc::as_ptr(&self.hap_to_seq) as usize)
    }
}

impl std::fmt::Display for HapRefGTRec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&to_vcf_rec(self))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ints::WrappedIntArray;
    use crate::vcf::MarkerParser;

    #[test]
    fn sequence_coded_access() {
        let marker = Marker::instance(
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1",
            &MarkerParser::new(true, true, true, true),
        );
        let samples = Samples::new(&["s0".into(), "s1".into()], &[true, true]);
        let hap_to_seq = WrappedIntArray::from_slice(&[0, 1, 0, 1]);
        let seq_to_allele = WrappedIntArray::from_slice(&[0, 1]);
        let r = HapRefGTRec::new(marker, samples, &hap_to_seq, &seq_to_allele);
        assert_eq!(r.size(), 4);
        assert_eq!(
            (0..4).map(|h| r.get(h)).collect::<Vec<_>>(),
            vec![0, 1, 0, 1]
        );
        assert_eq!(r.n_maps(), 2);
        let m0 = r.map(0);
        let m1 = r.map(1);
        assert_eq!(
            (0..4).map(|h| m0.get(h)).collect::<Vec<_>>(),
            vec![0, 1, 0, 1]
        );
        assert_eq!((0..2).map(|s| m1.get(s)).collect::<Vec<_>>(), vec![0, 1]);
        assert_eq!(r.allele_counts(), vec![2, 2]);
        assert_eq!(r.major_allele(), 0);
        let h2a = r.hap_to_allele();
        for h in 0..4 {
            assert_eq!(h2a.get(h), r.get(h));
        }
        assert_eq!(r.allele_to_haps()[1], Some(vec![1, 3]));
    }
}
