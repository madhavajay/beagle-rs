//! Port of `imp/HaplotypeCoder.java` — indexes the distinct allele sequences carried by the
//! phased reference and target haplotypes over a marker interval, returning, per interval, an
//! `IndexArray` mapping each haplotype (reference haps first, then target haps) to its sequence.

use std::rc::Rc;

use crate::ints::{IndexArray, IntArray};
use crate::vcf::{RefGT, GT};

/// Port of `imp/HaplotypeCoder.java`.
pub struct HaplotypeCoder {
    n_ref_haps: i32,
    n_haps: i32,
    ref_gt: RefGT,
    targ: Rc<dyn GT>,
}

/// The reference portion of a combined hap→sequence map.
enum RefPart {
    /// `codedRef[index]`.
    Direct(Vec<i32>),
    /// `seq1ToSeq2[refBasicIndexArray.get(index)]`.
    ViaSeq {
        basic: Box<dyn IntArray>,
        seq1_to_seq2: Vec<i32>,
    },
}

/// The lazily-evaluated combined hap→sequence map (Java's anonymous `IntArray`).
struct Combined {
    n_ref_haps: i32,
    n_haps: i32,
    ref_part: RefPart,
    targ_hap_to_seq: Vec<i32>,
}

impl IntArray for Combined {
    fn size(&self) -> i32 {
        self.n_haps
    }

    fn get(&self, index: i32) -> i32 {
        if index < self.n_ref_haps {
            match &self.ref_part {
                RefPart::Direct(coded_ref) => coded_ref[index as usize],
                RefPart::ViaSeq {
                    basic,
                    seq1_to_seq2,
                } => seq1_to_seq2[basic.get(index) as usize],
            }
        } else {
            self.targ_hap_to_seq[(index - self.n_ref_haps) as usize]
        }
    }
}

impl HaplotypeCoder {
    /// `new HaplotypeCoder(RefGT restrictRefGT, GT phasedTarg)`.
    pub fn new(restrict_ref_gt: RefGT, phased_targ: Rc<dyn GT>) -> Self {
        assert!(
            restrict_ref_gt.markers() == phased_targ.markers(),
            "inconsistent markers"
        );
        let n_ref_haps = restrict_ref_gt.n_haps();
        let n_haps = n_ref_haps + phased_targ.n_haps();
        HaplotypeCoder {
            n_ref_haps,
            n_haps,
            ref_gt: restrict_ref_gt,
            targ: phased_targ,
        }
    }

    /// `refHapPairs()`.
    pub fn ref_hap_pairs(&self) -> &RefGT {
        &self.ref_gt
    }

    /// `targHapPairs()`.
    pub fn targ_hap_pairs(&self) -> &Rc<dyn GT> {
        &self.targ
    }

    /// `run(int start, int end)` — hap→sequence-index map over markers `[start, end)`.
    pub fn run(&self, start: i32, end: i32) -> IndexArray {
        assert!(start < end, "start >= end");
        if is_hap_coded(&self.ref_gt, start, end) {
            self.code_seq_coded_ref(start, end)
        } else {
            self.code_seq(start, end)
        }
    }

    /// `codeTarg` — returns `(hapToSeq, seqCnt, seqMap)`.
    fn code_targ(&self, start: i32, len: usize) -> (Vec<i32>, i32, Vec<Vec<i32>>) {
        let n_targ_haps = self.targ.n_haps();
        let mut hap_to_seq = vec![1i32; n_targ_haps as usize];
        let mut seq_map: Vec<Vec<i32>> = vec![Vec::new(); len];
        let mut seq_cnt = 2;
        for (j, slot) in seq_map.iter_mut().enumerate() {
            let m = start + j as i32;
            let n_alleles = self.ref_gt.marker(m).n_alleles();
            *slot = vec![0i32; (seq_cnt * n_alleles) as usize];
            seq_cnt = 1;
            for h in 0..n_targ_haps {
                let index = n_alleles * hap_to_seq[h as usize] + self.targ.allele(m, h);
                if slot[index as usize] == 0 {
                    slot[index as usize] = seq_cnt;
                    seq_cnt += 1;
                }
                hap_to_seq[h as usize] = slot[index as usize];
            }
        }
        (hap_to_seq, seq_cnt, seq_map)
    }

    // `s`/`h` also index method calls (`seq_to_allele.get(s)`, `ref_gt.allele(m, h)`) that
    // clippy can't see; the explicit index is required.
    #[allow(clippy::needless_range_loop)]
    fn code_seq_coded_ref(&self, start: i32, end: i32) -> IndexArray {
        let (targ_hap_to_seq, seq_cnt, seq_map) = self.code_targ(start, (end - start) as usize);
        let start_rec = self.ref_gt.get(start);
        debug_assert_eq!(start_rec.n_maps(), 2);
        let mut seq_to_allele = start_rec.map(1);
        let mut seq1_to_seq2 = vec![1i32; seq_to_allele.size() as usize];
        for (j, slot) in seq_map.iter().enumerate() {
            let m = start + j as i32;
            let n_alleles = self.ref_gt.marker(m).n_alleles();
            let rec = self.ref_gt.get(m);
            debug_assert_eq!(rec.n_maps(), 2);
            seq_to_allele = rec.map(1);
            for s in 0..seq1_to_seq2.len() {
                if seq1_to_seq2[s] > 0 {
                    let index = seq1_to_seq2[s] * n_alleles + seq_to_allele.get(s as i32);
                    seq1_to_seq2[s] = slot[index as usize];
                }
            }
        }
        let basic = self.ref_gt.get(start).map(0);
        let combined = Combined {
            n_ref_haps: self.n_ref_haps,
            n_haps: self.n_haps,
            ref_part: RefPart::ViaSeq {
                basic,
                seq1_to_seq2,
            },
            targ_hap_to_seq,
        };
        IndexArray::from_int_array(Box::new(combined), seq_cnt)
    }

    #[allow(clippy::needless_range_loop)]
    fn code_seq(&self, start: i32, end: i32) -> IndexArray {
        let (targ_hap_to_seq, seq_cnt, seq_map) = self.code_targ(start, (end - start) as usize);
        let mut coded_ref_hap = vec![1i32; self.n_ref_haps as usize];
        for (j, slot) in seq_map.iter().enumerate() {
            let m = start + j as i32;
            let n_alleles = self.ref_gt.marker(m).n_alleles();
            for h in 0..coded_ref_hap.len() {
                if coded_ref_hap[h] > 0 {
                    let index = coded_ref_hap[h] * n_alleles + self.ref_gt.allele(m, h as i32);
                    coded_ref_hap[h] = slot[index as usize];
                }
            }
        }
        let combined = Combined {
            n_ref_haps: self.n_ref_haps,
            n_haps: self.n_haps,
            ref_part: RefPart::Direct(coded_ref_hap),
            targ_hap_to_seq,
        };
        IndexArray::from_int_array(Box::new(combined), seq_cnt)
    }
}

fn is_hap_coded(ref_gt: &RefGT, start: i32, end: i32) -> bool {
    let start_rec = ref_gt.get(start);
    if start_rec.is_allele_coded() {
        false
    } else {
        // Java compares the `map(0)` IntArray object; `seq_block_key` is that identity.
        let key = start_rec.seq_block_key();
        for m in (start + 1)..end {
            let rec = ref_gt.get(m);
            if rec.is_allele_coded() || rec.seq_block_key() != key {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::{
        allele_ref_gt_rec_from_parser, BasicGT, BasicGTRec, GTRec, MarkerParser, RefGTRec, Samples,
        VcfHeader, VcfRecGTParser, HEADER_PREFIX,
    };

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
    fn codes_ref_and_targ_haplotype_sequences() {
        // reference: 3 samples (6 haps); target: 2 samples (4 haps); same 3 markers.
        let hr = header(3);
        let ht = header(2);
        let mp = MarkerParser::new(true, true, true, true);
        let ref_lines = [
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1",
            "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|1\t1|0\t0|0",
            "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t1|1\t0|0\t0|1",
        ];
        let ref_recs: Vec<Rc<dyn RefGTRec>> = ref_lines
            .iter()
            .map(|l| {
                Rc::from(allele_ref_gt_rec_from_parser(&VcfRecGTParser::new(
                    &hr, l, &mp,
                )))
            })
            .collect();
        let ref_gt = RefGT::from_recs(ref_recs);

        let targ_lines = [
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|1",
            "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|0\t1|0",
            "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0",
        ];
        let targ_recs: Vec<Rc<dyn GTRec>> = targ_lines
            .iter()
            .map(|l| {
                Rc::new(BasicGTRec::from_parser(&VcfRecGTParser::new(&ht, l, &mp))) as Rc<dyn GTRec>
            })
            .collect();
        let targ: Rc<dyn GT> = Rc::new(BasicGT::new(targ_recs));

        let coder = HaplotypeCoder::new(ref_gt, targ);
        let ia = coder.run(0, 3);
        let arr = ia.int_array();
        assert_eq!(arr.size(), 10); // 6 ref haps + 4 targ haps
        let vs = ia.value_size();
        for h in 0..10 {
            let seq = arr.get(h);
            assert!(seq >= 0 && seq < vs, "hap {h} seq {seq} vs {vs}");
        }
        // identical reference haplotypes map to the same sequence index
        let _ = Samples::new(&["x".into()], &[true]); // touch import
    }
}
