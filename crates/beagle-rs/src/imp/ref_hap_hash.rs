//! Port of `imp/RefHapHash.java` — for the sublist of reference haplotypes that carry stored
//! HMM-state probability at a target marker, stores a hash code (from the rare-allele sequence
//! over a reference-marker interval) and the per-haplotype ALT alleles, for fast IBS matching.

use crate::ints::IntList;
use crate::jdk::Random;
use crate::vcf::{RefGT, GT};

use super::StateProbs;

/// Port of `imp/RefHapHash.java`.
pub struct RefHapHash {
    targ_marker: i32,
    ref_gt: RefGT,
    i2hap: Vec<i32>,
    alt_alleles: Vec<IntList>, // (marker offset, ALT allele) pairs per sublist haplotype
    i2hash: Vec<i32>,
    start: i32,
    end: i32,
}

fn i2hap(state_probs: &[Box<dyn StateProbs>], targ_marker: i32) -> Vec<i32> {
    let mut list: Vec<i32> = Vec::with_capacity(10 * state_probs.len());
    for sp in state_probs {
        let m = sp.n_states(targ_marker);
        for k in 0..m {
            list.push(sp.ref_hap(targ_marker, k));
        }
    }
    list.sort_unstable();
    list.dedup();
    list
}

impl RefHapHash {
    /// `new RefHapHash(AtomicReferenceArray<StateProbs>, int targMarker, RefGT, int start, int end)`.
    pub fn new(
        state_probs: &[Box<dyn StateProbs>],
        targ_marker: i32,
        ref_hap_pairs: RefGT,
        start: i32,
        end: i32,
    ) -> Self {
        assert!(start >= 0 && start < end, "{start}");
        assert!(end <= ref_hap_pairs.n_markers(), "{end}");
        let i2hap = i2hap(state_probs, targ_marker);
        let i2hash = vec![0i32; i2hap.len()];
        let alt_alleles: Vec<IntList> = (0..i2hap.len())
            .map(|_| IntList::with_capacity(6))
            .collect();
        let mut h = RefHapHash {
            targ_marker,
            ref_gt: ref_hap_pairs,
            i2hap,
            alt_alleles,
            i2hash,
            start,
            end,
        };
        h.set_hash_and_alt_alleles();
        h
    }

    fn set_hash_and_alt_alleles(&mut self) {
        let mut rand = Random::new(self.start as i64);
        for m in self.start..self.end {
            let rec = self.ref_gt.get(m);
            let marker_offset = m - self.start;
            if rec.is_allele_coded() && rec.major_allele() == 0 {
                self.low_alt_freq_update(rec.as_ref(), marker_offset, &mut rand);
            } else {
                self.standard_update(rec.as_ref(), marker_offset, &mut rand);
            }
        }
    }

    fn low_alt_freq_update(
        &mut self,
        rec: &dyn crate::vcf::RefGTRec,
        marker_offset: i32,
        rand: &mut Random,
    ) {
        let n_alleles = rec.marker().n_alleles();
        debug_assert_eq!(rec.major_allele(), 0);
        let n_haps = self.i2hap.len() as i32;
        for al in 1..n_alleles {
            let hash = rand.next_int();
            let n_copies = rec.allele_count(al);
            if n_haps < n_copies {
                for i in 0..self.i2hap.len() {
                    if rec.is_carrier(al, self.i2hap[i]) {
                        self.i2hash[i] = self.i2hash[i].wrapping_add(hash);
                        self.alt_alleles[i].add(marker_offset);
                        self.alt_alleles[i].add(al);
                    }
                }
            } else {
                for c in 0..n_copies {
                    let hap = rec.hap_index(al, c);
                    if let Ok(i) = self.i2hap.binary_search(&hap) {
                        self.i2hash[i] = self.i2hash[i].wrapping_add(hash);
                        self.alt_alleles[i].add(marker_offset);
                        self.alt_alleles[i].add(al);
                    }
                }
            }
        }
    }

    fn standard_update(
        &mut self,
        rec: &dyn crate::vcf::RefGTRec,
        marker_offset: i32,
        rand: &mut Random,
    ) {
        let n_alleles = rec.marker().n_alleles();
        let mut allele_hash = vec![0i32; n_alleles as usize];
        for slot in allele_hash.iter_mut().skip(1) {
            *slot = rand.next_int();
        }
        for i in 0..self.i2hap.len() {
            let allele = rec.get(self.i2hap[i]);
            if allele != 0 {
                self.i2hash[i] = self.i2hash[i].wrapping_add(allele_hash[allele as usize]);
                self.alt_alleles[i].add(marker_offset);
                self.alt_alleles[i].add(allele);
            }
        }
    }

    /// `targMarker()`.
    pub fn targ_marker(&self) -> i32 {
        self.targ_marker
    }
    /// `refGT()`.
    pub fn ref_gt(&self) -> &RefGT {
        &self.ref_gt
    }
    /// `start()`.
    pub fn start(&self) -> i32 {
        self.start
    }
    /// `end()`.
    pub fn end(&self) -> i32 {
        self.end
    }
    /// `nHaps()` — size of the reference-haplotype sublist.
    pub fn n_haps(&self) -> i32 {
        self.i2hap.len() as i32
    }
    /// `hap(int index)`.
    pub fn hap(&self, index: i32) -> i32 {
        self.i2hap[index as usize]
    }
    /// `hash(int index)`.
    pub fn hash(&self, index: i32) -> i32 {
        self.i2hash[index as usize]
    }

    /// `hap2Index(int hap)` — Java `Arrays.binarySearch` semantics (`-(insertionPoint) - 1`
    /// when absent).
    pub fn hap2_index(&self, hap: i32) -> i32 {
        match self.i2hap.binary_search(&hap) {
            Ok(i) => i as i32,
            Err(ip) => -(ip as i32) - 1,
        }
    }

    /// `setAlleles(int index, int[] alleles)` — writes the sublist haplotype's allele sequence
    /// over `[start, end)` into `alleles[0..end-start]` (0 = major allele).
    pub fn set_alleles(&self, index: i32, alleles: &mut [i32]) {
        let len = (self.end - self.start) as usize;
        alleles[..len].fill(0);
        let il = &self.alt_alleles[index as usize];
        let mut j = 0;
        while j < il.size() {
            alleles[il.get(j) as usize] = il.get(j + 1);
            j += 2;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::imp::StateProbsFactory;
    use crate::vcf::{
        allele_ref_gt_rec_from_parser, MarkerParser, RefGTRec, VcfHeader, VcfRecGTParser, GT,
        HEADER_PREFIX,
    };
    use std::rc::Rc;

    fn header(n: usize) -> VcfHeader {
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..n {
            hdr.push_str(&format!("\tS{s}"));
        }
        VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &vec![true; n],
        )
    }

    #[test]
    fn hashes_and_recovers_allele_sequences() {
        // 4 ref samples (8 haps), 3 markers
        let h = header(4);
        let mp = MarkerParser::new(true, true, true, true);
        let lines = [
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t0|0\t0|0",
            "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|1\t1|0\t0|0\t0|1",
            "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|0\t1|1\t0|0",
        ];
        let recs: Vec<Rc<dyn RefGTRec>> = lines
            .iter()
            .map(|l| {
                Rc::from(allele_ref_gt_rec_from_parser(&VcfRecGTParser::new(
                    &h, l, &mp,
                )))
            })
            .collect();
        let ref_gt = RefGT::from_recs(recs);

        // a StateProbs that references ref haps 0,1,2,3 at target marker 0
        let factory = StateProbsFactory::new(1);
        let hap_indices = vec![vec![0, 1, 2, 3]];
        let probs = vec![vec![0.4f32, 0.3, 0.2, 0.1]];
        let sp: Box<dyn StateProbs> = Box::new(factory.state_probs(0, 4, &hap_indices, &probs));
        let state_probs = vec![sp];

        let rhh = RefHapHash::new(&state_probs, 0, ref_gt, 0, 3);
        assert_eq!(rhh.start(), 0);
        assert_eq!(rhh.end(), 3);
        assert!(rhh.n_haps() >= 1);

        // hap2Index round-trips for present haps and is negative for absent
        for i in 0..rhh.n_haps() {
            let hap = rhh.hap(i);
            assert_eq!(rhh.hap2_index(hap), i);
        }
        assert!(rhh.hap2_index(99) < 0);

        // setAlleles recovers the ALT alleles for each sublist haplotype
        let mut alleles = vec![0i32; 3];
        for i in 0..rhh.n_haps() {
            rhh.set_alleles(i, &mut alleles);
            let hap = rhh.hap(i);
            for (m, &al) in alleles.iter().enumerate() {
                assert_eq!(al, rhh.ref_gt().allele(m as i32, hap), "i={i} m={m}");
            }
        }
    }
}
