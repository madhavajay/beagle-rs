//! Port of `phase/Stage2Baum.java` — runs the haploid Li & Stephens forward-backward HMM at
//! the high-frequency (stage-1) markers and uses the resulting state probabilities to impute
//! missing genotypes and resolve heterozygote phase at the low-frequency markers.

use std::rc::Rc;

use crate::ints::IntArray;
use crate::jdk::Random;
use crate::vcf::{RefGT, GT};

use super::{FixedPhaseData, HmmStateProbs, LowFreqPhaseIbs, Stage2Haps};

/// Port of `phase/Stage2Baum.java`.
pub struct Stage2Baum<'a, 'b> {
    fpd: Rc<FixedPhaseData>,
    state_probs: HmmStateProbs<'a>,
    n_states: [i32; 2],
    states: Vec<Vec<Vec<i32>>>, // [2][nStage1Markers][maxStates]
    probs: Vec<Vec<Vec<f32>>>,  // [2][nStage1Markers][maxStates]
    unph_targ_gt: Rc<dyn GT>,
    ref_gt: Option<RefGT>,
    n_targ_haps: i32,
    n_stage1_markers: i32,
    stage2_haps: &'b mut Stage2Haps,
    stage1_to2: crate::ints::WrappedIntArray,
    seed: i64,
    rand: Random,
}

fn max_index(fa: &[f32]) -> i32 {
    let mut max_index = 0usize;
    for j in 1..fa.len() {
        if fa[j] > fa[max_index] {
            max_index = j;
        }
    }
    max_index as i32
}

impl<'a, 'b> Stage2Baum<'a, 'b> {
    /// `new Stage2Baum(LowFreqPhaseIbs phaseIbs, Stage2Haps stage2Haps)`.
    pub fn new(phase_ibs: &'a LowFreqPhaseIbs<'a>, stage2_haps: &'b mut Stage2Haps) -> Self {
        let phase_data = phase_ibs.phase_data();
        let fpd = phase_data.est_phase().fpd_rc().clone();
        let n_stage1_markers = fpd.stage1_targ_gt().n_markers();
        let state_probs = HmmStateProbs::new(phase_ibs);
        let max_states = state_probs.max_states() as usize;
        let states = vec![vec![vec![0i32; max_states]; n_stage1_markers as usize]; 2];
        let probs = vec![vec![vec![0.0f32; max_states]; n_stage1_markers as usize]; 2];
        let unph_targ_gt = fpd.targ_gt().clone();
        let ref_gt = fpd.restricted_ref_gt().cloned();
        let n_targ_haps = fpd.targ_gt().n_haps();
        let stage1_to2 = fpd.stage1_to2().clone();
        let seed = phase_data.seed();
        Stage2Baum {
            fpd,
            state_probs,
            n_states: [0, 0],
            states,
            probs,
            unph_targ_gt,
            ref_gt,
            n_targ_haps,
            n_stage1_markers,
            stage2_haps,
            stage1_to2,
            seed,
            rand: Random::new(seed),
        }
    }

    /// `nTargSamples()`.
    pub fn n_targ_samples(&self) -> i32 {
        self.fpd.targ_gt().n_samples()
    }

    /// `phase(int targSample)`.
    pub fn phase(&mut self, targ_sample: i32) {
        self.rand.set_seed(self.seed + targ_sample as i64);
        let h1 = targ_sample << 1;
        let h2 = h1 | 0b1;
        self.n_states[0] = self
            .state_probs
            .run(h1, &mut self.states[0], &mut self.probs[0]);
        self.n_states[1] = self
            .state_probs
            .run(h2, &mut self.states[1], &mut self.probs[1]);
        let mut start = 0;
        for j in 0..self.n_stage1_markers {
            let end = self.stage1_to2.get(j);
            self.impute_interval(targ_sample, start, end);
            start = end + 1;
        }
        let n_markers = self.unph_targ_gt.n_markers();
        self.impute_interval(targ_sample, start, n_markers);
    }

    fn impute_interval(&mut self, sample: i32, start: i32, end: i32) {
        let hap1 = sample << 1;
        let hap2 = hap1 | 0b1;
        for m in start..end {
            let mut a1 = self.unph_targ_gt.allele(m, hap1);
            let mut a2 = self.unph_targ_gt.allele(m, hap2);
            if a1 >= 0 && a2 >= 0 {
                if a1 != a2 {
                    let al_probs1 = self.unscaled_al_probs(m, 0, a1, a2);
                    let al_probs2 = self.unscaled_al_probs(m, 1, a1, a2);
                    let p1 = al_probs1[a1 as usize] * al_probs2[a2 as usize];
                    let p2 = al_probs1[a2 as usize] * al_probs2[a1 as usize];
                    let switch_alleles = p1 < p2 || (p1 == p2 && self.rand.next_boolean());
                    if switch_alleles {
                        std::mem::swap(&mut a1, &mut a2);
                    }
                }
            } else {
                a1 = self.impute_allele(m, 0);
                a2 = self.impute_allele(m, 1);
            }
            self.stage2_haps.set_phased_gt(m, sample, a1, a2);
        }
    }

    fn unscaled_al_probs(&self, m: i32, hap_bit: i32, a1: i32, a2: i32) -> Vec<f32> {
        let n_alleles = self.unph_targ_gt.marker(m).n_alleles() as usize;
        let mut al_probs = vec![0.0f32; n_alleles];
        let rare1 = self.fpd.is_low_freq(m, a1);
        let rare2 = self.fpd.is_low_freq(m, a2);
        let mkr_a = self.fpd.prev_stage1_marker(m);
        let mkr_b = (mkr_a + 1).min(self.n_stage1_markers - 1);
        let hb = hap_bit as usize;
        let n = self.n_states[hb] as usize;
        for j in 0..n {
            let hap = self.states[hb][mkr_a as usize][j];
            let b1 = self.allele(m, hap);
            let b2 = self.allele(m, hap ^ 0b1);
            if b1 >= 0 && b2 >= 0 {
                let wt = self.fpd.prev_stage1_wt(m);
                let prob = wt * self.probs[hb][mkr_a as usize][j]
                    + (1.0 - wt) * self.probs[hb][mkr_b as usize][j];
                if b1 == b2 {
                    al_probs[b1 as usize] += prob;
                } else {
                    let match1 = rare1 && (a1 == b1 || a1 == b2);
                    let match2 = rare2 && (a2 == b1 || a2 == b2);
                    if match1 ^ match2 {
                        if match1 {
                            al_probs[a1 as usize] += prob;
                        } else {
                            al_probs[a2 as usize] += prob;
                        }
                    }
                }
            }
        }
        al_probs
    }

    fn impute_allele(&self, m: i32, hap_bit: i32) -> i32 {
        let n_alleles = self.unph_targ_gt.marker(m).n_alleles() as usize;
        let mut al_probs = vec![0.0f32; n_alleles];
        let mkr_a = self.fpd.prev_stage1_marker(m);
        let mkr_b = (mkr_a + 1).min(self.n_stage1_markers - 1);
        let hb = hap_bit as usize;
        let n = self.n_states[hb] as usize;
        for j in 0..n {
            let wt = self.fpd.prev_stage1_wt(m);
            let prob = wt * self.probs[hb][mkr_a as usize][j]
                + (1.0 - wt) * self.probs[hb][mkr_b as usize][j];
            let hap = self.states[hb][mkr_a as usize][j];
            let b1 = self.allele(m, hap);
            let b2 = self.allele(m, hap ^ 0b1);
            if b1 >= 0 && b2 >= 0 {
                if b1 == b2 || hap >= self.n_targ_haps {
                    al_probs[b1 as usize] += prob;
                } else {
                    let is_rare1 = self.fpd.is_low_freq(m, b1);
                    let is_rare2 = self.fpd.is_low_freq(m, b2);
                    // Java mixes 0.55/0.45/0.5 (double) with the float `prob`; reproduce the
                    // float←float+double widening/narrowing exactly.
                    if is_rare1 ^ is_rare2 {
                        if is_rare1 {
                            al_probs[b1 as usize] =
                                (al_probs[b1 as usize] as f64 + 0.55 * prob as f64) as f32;
                            al_probs[b2 as usize] =
                                (al_probs[b2 as usize] as f64 + 0.45 * prob as f64) as f32;
                        } else {
                            al_probs[b1 as usize] =
                                (al_probs[b1 as usize] as f64 + 0.45 * prob as f64) as f32;
                            al_probs[b2 as usize] =
                                (al_probs[b2 as usize] as f64 + 0.55 * prob as f64) as f32;
                        }
                    } else {
                        al_probs[b1 as usize] =
                            (al_probs[b1 as usize] as f64 + 0.5 * prob as f64) as f32;
                        al_probs[b2 as usize] =
                            (al_probs[b2 as usize] as f64 + 0.5 * prob as f64) as f32;
                    }
                }
            }
        }
        max_index(&al_probs)
    }

    fn allele(&self, marker: i32, hap: i32) -> i32 {
        if hap < self.n_targ_haps {
            self.unph_targ_gt.allele(marker, hap)
        } else {
            self.ref_gt
                .as_ref()
                .expect("ref")
                .allele(marker, hap - self.n_targ_haps)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::main_pkg::{Par, Pedigree};
    use crate::vcf::{
        BasicGT, BasicGTRec, GTRec, GeneticMap, MarkerIndices, MarkerParser, PositionMap,
        VcfHeader, VcfRecGTParser, Window, HEADER_PREFIX,
    };
    use std::io::Write;

    use super::super::PhaseData;

    fn fpd() -> FixedPhaseData {
        let gt = std::env::temp_dir().join("beagle_rs_s2b_gt.vcf");
        std::fs::File::create(&gt).unwrap().write_all(b"x").unwrap();
        let par = Par::new(&[
            format!("gt={}", gt.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
            "seed=99999".to_string(),
        ]);
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..6 {
            hdr.push_str(&format!("\tS{s}"));
        }
        let h = VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &[true; 6],
        );
        let mp = MarkerParser::new(true, true, true, true);
        let gts = [
            "0|0\t0|1\t1|1\t0|1\t0|0\t1|0",
            "0|1\t0|0\t1|0\t1|1\t0|1\t0|0",
            "1|1\t0|1\t0|0\t0|1\t1|0\t0|1",
            "0|0\t1|1\t0|1\t1|0\t0|0\t1|1",
            "0|1\t1|0\t1|1\t0|0\t0|1\t1|0",
            "1|0\t0|1\t0|0\t1|1\t1|0\t0|1",
        ];
        let recs: Vec<Rc<dyn GTRec>> = gts
            .iter()
            .enumerate()
            .map(|(i, g)| {
                let pos = 1_000_000 + i as i32 * 1_000_000;
                let line = format!("chr1\t{pos}\t.\tA\tC\t.\tPASS\t.\tGT\t{g}");
                Rc::new(BasicGTRec::from_parser(&VcfRecGTParser::new(
                    &h, &line, &mp,
                ))) as Rc<dyn GTRec>
            })
            .collect();
        let targ = BasicGT::new(recs);
        let n = targ.n_markers();
        let indices = MarkerIndices::from_counts(0, n, n);
        let gen_map: Rc<dyn GeneticMap> = Rc::new(PositionMap::new(1e-6));
        let w = Window::new(gen_map, 1, true, indices, None, targ);
        let ped = Pedigree::new(w.targ_gt().samples().clone(), None);
        FixedPhaseData::new(&par, &ped, &w, None)
    }

    #[test]
    fn phases_all_samples_into_stage2haps() {
        let pd = PhaseData::new(Rc::new(fpd()), 99999);
        let n_markers = pd.fpd().targ_gt().n_markers();
        let phase_ibs = LowFreqPhaseIbs::new(&pd);
        let mut stage2_haps = Stage2Haps::new(&pd);
        {
            let mut baum = Stage2Baum::new(&phase_ibs, &mut stage2_haps);
            assert_eq!(baum.n_targ_samples(), 6);
            for s in 0..baum.n_targ_samples() {
                baum.phase(s);
            }
        }
        let gt = stage2_haps.to_basic_gt(0, n_markers);
        assert_eq!(gt.n_markers(), n_markers);
        assert_eq!(gt.n_samples(), 6);
        // every genotype is now phased and non-missing
        for m in 0..n_markers {
            for hap in 0..12 {
                assert!(gt.allele(m, hap) >= 0, "m={m} hap={hap} unphased");
            }
        }
        let _ = &phase_ibs; // keep alive
    }
}
