//! Port of `bref/SeqCoder3.java` — compresses a run of allele-coded `RefGTRec`s by storing
//! the distinct allele sequences and each haplotype's sequence index (bref v3 seq-coding).

use std::rc::Rc;

use crate::ints::{packed_create, IntList};
use crate::vcf::{HapRefGTRec, RefGTRec, Samples};

/// `SeqCoder3.MAX_NALLELES`.
pub const MAX_NALLELES: i32 = 255;
/// `SeqCoder3.COMPRESS_FREQ_THRESHOLD`.
pub const COMPRESS_FREQ_THRESHOLD: f32 = 0.995;

const CHAR_MAX_VALUE: i32 = 65535; // Java `Character.MAX_VALUE`

/// Port of `bref/SeqCoder3.java`.
pub struct SeqCoder3 {
    samples: Samples,
    max_n_seq: i32,
    recs: Vec<Rc<dyn RefGTRec>>,
    hap2_seq: Vec<i32>,
    seq2_cnt: IntList,
    seq2_allele_seq_map: Vec<IntList>,
}

/// `SeqCoder3.defaultMaxNSeq(int nSamples)`.
pub fn default_max_n_seq(n_samples: i32) -> i32 {
    assert!(n_samples >= 1, "{n_samples}");
    if n_samples == 1 {
        3
    } else {
        let exponent = 2.0 * (n_samples as f64).log10() + 1.0;
        let max_n_seq = 2.0f64.powf(exponent).floor() as i64;
        if max_n_seq > CHAR_MAX_VALUE as i64 {
            CHAR_MAX_VALUE
        } else {
            max_n_seq as i32
        }
    }
}

impl SeqCoder3 {
    /// `new SeqCoder3(Samples samples)`.
    pub fn new(samples: Samples) -> Self {
        let max = default_max_n_seq(samples.size());
        SeqCoder3::with_max_n_seq(samples, max)
    }

    /// `new SeqCoder3(Samples samples, int maxNSeq)`.
    pub fn with_max_n_seq(samples: Samples, max_n_seq: i32) -> Self {
        assert!(max_n_seq >= 0 && max_n_seq < CHAR_MAX_VALUE, "{max_n_seq}");
        let n_haps = 2 * samples.size();
        let mut coder = SeqCoder3 {
            samples,
            max_n_seq,
            recs: Vec::with_capacity(100),
            hap2_seq: vec![0; n_haps as usize],
            seq2_cnt: IntList::with_capacity(3 * max_n_seq / 2 + 1),
            seq2_allele_seq_map: Vec::with_capacity((max_n_seq + 1) as usize),
        };
        coder.initialize();
        coder
    }

    /// `samples()`.
    pub fn samples(&self) -> &Samples {
        &self.samples
    }

    /// `nRecs()`.
    pub fn n_recs(&self) -> i32 {
        self.recs.len() as i32
    }

    /// `maxNSeq()`.
    pub fn max_n_seq(&self) -> i32 {
        self.max_n_seq
    }

    /// `add(RefGTRec rec)` — attempts to add `rec`; returns `true` if it was added.
    pub fn add(&mut self, rec: Rc<dyn RefGTRec>) -> bool {
        assert!(rec.samples() == &self.samples, "inconsistent samples");
        assert!(rec.is_allele_coded(), "record is not allele-coded");
        let success = self.set_allele_map(rec.as_ref());
        if success {
            let major_allele = rec.major_allele();
            let n_alleles = rec.marker().n_alleles();
            for a in 0..n_alleles {
                if a != major_allele {
                    let n_copies = rec.allele_count(a);
                    for c in 0..n_copies {
                        let h = rec.hap_index(a, c) as usize;
                        let old_seq = self.hap2_seq[h];
                        let list = &self.seq2_allele_seq_map[old_seq as usize];
                        let mut index = 0;
                        while index < list.size() && list.get(index) != a {
                            index += 2;
                        }
                        let new_seq = list.get(index + 1);
                        if new_seq != old_seq {
                            while new_seq >= self.seq2_cnt.size() {
                                self.seq2_cnt.add(0);
                            }
                            self.hap2_seq[h] = new_seq;
                            self.seq2_cnt.decrement_and_get(old_seq);
                            self.seq2_cnt.increment_and_get(new_seq);
                        }
                    }
                }
            }
            self.recs.push(rec);
        }
        debug_assert_eq!(
            self.seq2_cnt.size() as usize,
            self.seq2_allele_seq_map.len()
        );
        success
    }

    fn clear_seq2_allele_map(&mut self) {
        for list in &mut self.seq2_allele_seq_map {
            list.clear();
        }
    }

    fn set_allele_map(&mut self, rec: &dyn RefGTRec) -> bool {
        let n_start_seq = self.seq2_cnt.size();
        let mut seq2_non_major_cnt = vec![0i32; n_start_seq as usize];
        self.clear_seq2_allele_map();
        let n_alleles = rec.marker().n_alleles();
        let major_allele = rec.major_allele();
        for a in 0..n_alleles {
            if a != major_allele {
                let n_copies = rec.allele_count(a);
                for c in 0..n_copies {
                    let h = rec.hap_index(a, c) as usize;
                    let seq = self.hap2_seq[h] as usize;
                    seq2_non_major_cnt[seq] += 1;
                    if self.seq2_allele_seq_map[seq].is_empty() {
                        self.seq2_allele_seq_map[seq].add(a);
                        self.seq2_allele_seq_map[seq].add(seq as i32);
                    } else {
                        let mut index = 0;
                        while index < self.seq2_allele_seq_map[seq].size()
                            && self.seq2_allele_seq_map[seq].get(index) != a
                        {
                            index += 2;
                        }
                        if index == self.seq2_allele_seq_map[seq].size() {
                            let new_seq = self.seq2_allele_seq_map.len() as i32;
                            self.seq2_allele_seq_map[seq].add(a);
                            self.seq2_allele_seq_map[seq].add(new_seq);
                            self.seq2_allele_seq_map.push(IntList::with_capacity(4));
                        }
                    }
                }
            }
        }
        self.add_major_allele(&seq2_non_major_cnt, major_allele);
        if self.seq2_allele_seq_map.len() as i32 > self.max_n_seq {
            self.seq2_allele_seq_map.truncate(n_start_seq as usize);
            false
        } else {
            true
        }
    }

    fn add_major_allele(&mut self, seq2_non_major_cnt: &[i32], major_allele: i32) {
        for (seq, &non_major_cnt) in seq2_non_major_cnt.iter().enumerate() {
            if non_major_cnt < self.seq2_cnt.get(seq as i32) {
                if self.seq2_allele_seq_map[seq].is_empty() {
                    self.seq2_allele_seq_map[seq].add(major_allele);
                    self.seq2_allele_seq_map[seq].add(seq as i32);
                } else {
                    let g0 = self.seq2_allele_seq_map[seq].get(0);
                    debug_assert_eq!(self.seq2_allele_seq_map[seq].get(1), seq as i32);
                    let new_seq = self.seq2_allele_seq_map.len() as i32;
                    self.seq2_allele_seq_map[seq].add(g0);
                    self.seq2_allele_seq_map[seq].add(new_seq);
                    self.seq2_allele_seq_map[seq].set(0, major_allele);
                    self.seq2_allele_seq_map.push(IntList::with_capacity(4));
                }
            }
        }
    }

    /// `getCompressedList()` — returns the compressed `HapRefGTRec`s and resets the coder.
    pub fn get_compressed_list(&mut self) -> Vec<Box<dyn RefGTRec>> {
        if self.recs.is_empty() {
            return Vec::new();
        }
        let mut list: Vec<Box<dyn RefGTRec>> = Vec::with_capacity(self.recs.len());
        let seq2_hap = self.seq2_first_hap();
        // One shared `hap_to_seq` map for the whole group, so the records compare equal by
        // identity for bref3 block grouping (Java reuses one `IntArray` object per group).
        let hap2seq = std::rc::Rc::new(self.hap2_seq.clone());
        let recs = std::mem::take(&mut self.recs);
        for rec in &recs {
            let m = rec.marker().clone();
            let seq2allele = seq2_allele(rec.as_ref(), &seq2_hap);
            list.push(Box::new(HapRefGTRec::new_shared(
                m,
                self.samples.clone(),
                hap2seq.clone(),
                seq2allele.as_ref(),
            )));
        }
        self.initialize();
        list
    }

    fn seq2_first_hap(&self) -> Vec<i32> {
        let mut seq_to_first_hap = vec![-1i32; self.seq2_allele_seq_map.len()];
        for (h, &seq) in self.hap2_seq.iter().enumerate() {
            if seq_to_first_hap[seq as usize] == -1 {
                seq_to_first_hap[seq as usize] = h as i32;
            }
        }
        seq_to_first_hap
    }

    fn initialize(&mut self) {
        self.recs.clear();
        self.seq2_cnt.clear();
        self.seq2_allele_seq_map.clear();
        // initialize with the empty sequence (seq index 0)
        self.hap2_seq.iter_mut().for_each(|x| *x = 0);
        self.seq2_cnt.add(self.hap2_seq.len() as i32);
        self.seq2_allele_seq_map.push(IntList::with_capacity(4));
    }
}

fn seq2_allele(rec: &dyn RefGTRec, seq2_hap: &[i32]) -> Box<dyn crate::ints::IntArray> {
    let seq2_allele: Vec<i32> = seq2_hap.iter().map(|&h| rec.get(h)).collect();
    packed_create(&seq2_allele, rec.marker().n_alleles())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::{allele_ref_gt_rec_from_components, Marker, MarkerParser};

    fn marker(line: &str) -> Marker {
        Marker::instance(line, &MarkerParser::new(true, true, true, true))
    }

    fn rec(samples: &Samples, line: &str, alts: Vec<i32>) -> Rc<dyn RefGTRec> {
        // biallelic: allele 1 carried by `alts`, major allele 0 (None).
        Rc::from(allele_ref_gt_rec_from_components(
            marker(line),
            samples.clone(),
            vec![None, Some(alts)],
        ))
    }

    #[test]
    fn default_max_n_seq_values() {
        assert_eq!(default_max_n_seq(1), 3);
        // 2*log10(100)+1 = 5 -> 2^5 = 32
        assert_eq!(default_max_n_seq(100), 32);
    }

    #[test]
    fn compress_round_trip_reproduces_genotypes() {
        let samples = Samples::new(
            &["s0".into(), "s1".into(), "s2".into(), "s3".into()],
            &[true; 4],
        );
        let r0 = rec(
            &samples,
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1",
            vec![1, 3],
        );
        let r1 = rec(
            &samples,
            "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|1",
            vec![0, 1],
        );
        let n_haps = samples.size() << 1;

        // expected per-hap alleles from the originals
        let exp0: Vec<i32> = (0..n_haps).map(|h| r0.get(h)).collect();
        let exp1: Vec<i32> = (0..n_haps).map(|h| r1.get(h)).collect();

        let mut coder = SeqCoder3::with_max_n_seq(samples, 100);
        assert!(coder.add(r0));
        assert!(coder.add(r1));
        assert_eq!(coder.n_recs(), 2);

        let compressed = coder.get_compressed_list();
        assert_eq!(compressed.len(), 2);
        assert_eq!(coder.n_recs(), 0); // reset after flush
        let got0: Vec<i32> = (0..n_haps).map(|h| compressed[0].get(h)).collect();
        let got1: Vec<i32> = (0..n_haps).map(|h| compressed[1].get(h)).collect();
        assert_eq!(got0, exp0);
        assert_eq!(got1, exp1);
        // seq-coded records are phased reference data
        assert!(compressed[0].is_phased());
    }
}
