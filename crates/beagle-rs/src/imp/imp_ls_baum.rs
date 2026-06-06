//! Port of `imp/ImpLSBaum.java` — the Li & Stephens forward-backward HMM that computes, for a
//! target haplotype, the posterior HMM state probabilities at the genotyped marker clusters
//! using the `ImpStates` IBS-matched reference haplotypes.

use super::{BasicStateProbs, ImpData, ImpIbs, ImpStates, StateProbsFactory};

/// Port of `imp/ImpLSBaum.java`.
pub struct ImpLSBaum<'a> {
    imp_data: &'a ImpData,
    states: ImpStates<'a>,
    n_markers: i32,
    hap_indices: Vec<Vec<i32>>,
    alleles_match: Vec<Vec<bool>>,
    fwd_val: Vec<Vec<f32>>,
    bwd_val: Vec<f32>,
    al_probs_factory: StateProbsFactory,
    targ_allele: Vec<i32>,
}

impl<'a> ImpLSBaum<'a> {
    /// `new ImpLSBaum(ImpData impData, ImpIbs ibsHaps)`.
    pub fn new(imp_data: &'a ImpData, ibs_haps: &'a ImpIbs<'a>) -> Self {
        let states = ImpStates::new(ibs_haps);
        let n_markers = imp_data.n_clusters();
        let max_states = imp_data.par().imp_states() as usize;
        let nm = n_markers as usize;
        ImpLSBaum {
            imp_data,
            states,
            n_markers,
            hap_indices: vec![vec![0i32; max_states]; nm],
            alleles_match: vec![vec![false; max_states]; nm],
            fwd_val: vec![vec![0f32; max_states]; nm],
            bwd_val: vec![0f32; max_states],
            al_probs_factory: StateProbsFactory::new(n_markers),
            targ_allele: vec![0i32; imp_data.n_clusters() as usize],
        }
    }

    /// `impData()`.
    pub fn imp_data(&self) -> &ImpData {
        self.imp_data
    }

    /// `impute(int targHap)` — posterior HMM state probabilities at the genotyped markers.
    pub fn impute(&mut self, targ_hap: i32) -> BasicStateProbs {
        let last_marker = self.imp_data.n_clusters() - 1;
        let n_states =
            self.states
                .ibs_states(targ_hap, &mut self.hap_indices, &mut self.alleles_match);
        self.set_fwd_values(targ_hap, n_states);
        self.bwd_val[..n_states as usize].fill(1.0 / n_states as f32);
        let mut last_sum = 1.0f32;
        for m in (0..=last_marker).rev() {
            last_sum = self.set_bwd_value(m, n_states, last_sum);
        }
        self.al_probs_factory
            .state_probs(targ_hap, n_states, &self.hap_indices, &self.fwd_val)
    }

    fn set_fwd_values(&mut self, targ_hap: i32, n_haps: i32) {
        let n_ref_haps = self.imp_data.n_ref_haps();
        let mut last_sum = 1.0f32;
        let len = self.fwd_val.len();
        let n = n_haps as usize;
        for m in 0..len {
            let p_recomb = self.imp_data.p_recomb(m as i32);
            let p_err = self.imp_data.err_prob(m as i32);
            let p_no_err = 1.0 - p_err;
            let shift = p_recomb / n_haps as f32;
            let scale = (1.0 - p_recomb) / last_sum;
            let mut sum = 0.0f32;
            self.targ_allele[m] = self.imp_data.allele(m as i32, n_ref_haps + targ_hap);
            if m == 0 {
                let am = &self.alleles_match[0];
                let cur = &mut self.fwd_val[0];
                for j in 0..n {
                    let em = if am[j] { p_no_err } else { p_err };
                    cur[j] = em;
                    sum += cur[j];
                }
            } else {
                let (head, tail) = self.fwd_val.split_at_mut(m);
                let prev = &head[m - 1];
                let cur = &mut tail[0];
                let am = &self.alleles_match[m];
                for j in 0..n {
                    let em = if am[j] { p_no_err } else { p_err };
                    cur[j] = em * (scale * prev[j] + shift);
                    sum += cur[j];
                }
            }
            last_sum = sum;
        }
    }

    fn set_bwd_value(&mut self, m: i32, n_states: i32, last_sum: f32) -> f32 {
        let m_p1 = m + 1;
        let p_recomb = if m_p1 < self.n_markers {
            self.imp_data.p_recomb(m_p1)
        } else {
            0.0
        };
        let p_err = self.imp_data.err_prob(m);
        let p_no_err = 1.0 - p_err;
        let scale = (1.0 - p_recomb) / last_sum;
        let shift = p_recomb / n_states as f32;
        let mut bwd_val_sum = 0f32;
        let mut state_sum = 0f32;
        let ns = n_states as usize;
        let am = &self.alleles_match[m as usize];
        let fwd_m = &mut self.fwd_val[m as usize];
        for j in 0..ns {
            self.bwd_val[j] = scale * self.bwd_val[j] + shift; // finish bwd value
            fwd_m[j] *= self.bwd_val[j]; // store state probability in fwd_val[m]
            state_sum += fwd_m[j];
            let em = if am[j] { p_no_err } else { p_err };
            self.bwd_val[j] *= em;
            bwd_val_sum += self.bwd_val[j];
        }
        for v in fwd_m.iter_mut().take(ns) {
            *v /= state_sum; // normalize state probabilities
        }
        bwd_val_sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::imp::StateProbs;
    use crate::main_pkg::Par;
    use crate::vcf::{
        allele_ref_gt_rec_from_parser, BasicGT, BasicGTRec, GTRec, GeneticMap, MarkerIndices,
        MarkerParser, PositionMap, RefGT, RefGTRec, VcfHeader, VcfRecGTParser, Window, GT,
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

    fn imp_data() -> ImpData {
        let gt = std::env::temp_dir().join("beagle_rs_implsbaum_gt.vcf");
        let rf = std::env::temp_dir().join("beagle_rs_implsbaum_ref.vcf");
        std::fs::write(&gt, b"x").unwrap();
        std::fs::write(&rf, b"x").unwrap();
        let par = Par::new(&[
            format!("gt={}", gt.display()),
            format!("ref={}", rf.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
        ]);
        let hr = header(6);
        let ht = header(3);
        let mp = MarkerParser::new(true, true, true, true);
        let ref_lines = [
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\t0|1\t1|0\t0|0",
            "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|1\t1|0\t0|0\t1|1\t0|1\t1|0",
            "chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t1|1\t0|0\t0|1\t0|0\t1|0\t0|1",
            "chr1\t400\t.\tC\tA\t.\tPASS\t.\tGT\t0|1\t0|1\t1|0\t1|1\t0|0\t1|0",
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
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|1\t0|0",
            "chr1\t400\t.\tC\tA\t.\tPASS\t.\tGT\t0|0\t1|0\t0|1",
        ];
        let targ_recs: Vec<Rc<dyn GTRec>> = targ_lines
            .iter()
            .map(|l| {
                Rc::new(BasicGTRec::from_parser(&VcfRecGTParser::new(&ht, l, &mp))) as Rc<dyn GTRec>
            })
            .collect();
        let targ = BasicGT::new(targ_recs);
        let indices = MarkerIndices::from_in_targ(&[true, false, false, true], 0, 4);
        let gen_map: Rc<dyn GeneticMap> = Rc::new(PositionMap::new(1e-6));
        let window = Window::new(gen_map.clone(), 1, true, indices, Some(ref_gt), targ);
        let phased: Rc<dyn GT> = Rc::new(window.targ_gt().clone());
        ImpData::new(&par, &window, phased, gen_map.as_ref())
    }

    #[test]
    fn imputes_state_probabilities() {
        let imp = imp_data();
        let n_clusters = imp.n_clusters();
        let ibs = ImpIbs::new(&imp);
        let mut baum = ImpLSBaum::new(&imp, &ibs);
        let sp = baum.impute(0);
        assert_eq!(sp.targ_hap(), 0);
        assert_eq!(sp.n_targ_markers(), n_clusters);
        for m in 0..n_clusters {
            for k in 0..sp.n_states(m) {
                let p = sp.probs(m, k);
                assert!((0.0..=1.0001).contains(&p), "m={m} k={k} p={p}");
                assert!(sp.ref_hap(m, k) >= 0);
            }
        }
    }
}
