//! Port of `phase/PhaseBaum2.java` — the stage-1 Li & Stephens forward-backward HMM that
//! re-phases unphased heterozygotes, imputes missing genotypes, and decides masked-het phase
//! for one target sample per call.
//!
//! Java mutates the single shared `SamplePhase` object obtained from `EstPhase`; the Rust port
//! works on `MarkerCluster`'s owned clone and stores it back (see [[estphase-get-set-aliasing]]
//! for why the masked clone is written before constructing the `MarkerCluster`).

use crate::ints::IntArray;
use crate::vcf::Markers;

use super::{
    BasicPhaseStates, ClustType, HmmUpdater, MarkerCluster, PbwtPhaseIbs, PhaseBaum, SamplePhase,
    SwapRate,
};

/// Port of `phase/PhaseBaum2.java`.
pub struct PhaseBaum2<'a> {
    phase_ibs: &'a PbwtPhaseIbs<'a>,
    burnin: bool,
    lr_threshold: f32,
    mask_trailing_hets: bool,
    markers: Markers,
    ref_alleles: Vec<Vec<i32>>, // ref panel alleles for each missing genotype
    mismatches: Vec<Vec<Vec<u8>>>, // [3][nMarker][nStates]
    p_mismatch: f32,
    em_probs: [f32; 2],
    max_states: i32,
    states: BasicPhaseStates<'a>,
    n_states: i32,
    fwd: [Vec<f32>; 3],
    bwd: [Vec<f32>; 3],
    fwd_sums: [f32; 3],
    bwd_miss1: Vec<Vec<f32>>,
    bwd_miss2: Vec<Vec<f32>>,
    bwd_het1: Vec<Vec<f32>>,
    bwd_het2: Vec<Vec<f32>>,
    swap_haps: bool,
    n_swaps: i32,
}

impl<'a> PhaseBaum2<'a> {
    /// `new PhaseBaum2(PbwtPhaseIbs phaseIbs)`.
    pub fn new(phase_ibs: &'a PbwtPhaseIbs<'a>) -> Self {
        let phase_data = phase_ibs.phase_data();
        let burnin = phase_data.it() < phase_data.fpd().par().burnin();
        let lr_threshold = phase_data.lr_threshold();
        let lr_mask_threshold = 50f32;
        let mask_trailing_hets = lr_threshold < lr_mask_threshold;
        let markers = phase_data.fpd().stage1_targ_gt().markers().clone();
        let n_markers = markers.size(); // only used to size `mismatches`
        let max_states = phase_data.fpd().par().phase_states();
        let states = BasicPhaseStates::new(phase_ibs, max_states);
        let p_mismatch = phase_data.p_mismatch();
        let ms = max_states as usize;
        PhaseBaum2 {
            phase_ibs,
            burnin,
            lr_threshold,
            mask_trailing_hets,
            markers,
            ref_alleles: Vec::new(),
            mismatches: vec![vec![vec![0u8; ms]; n_markers as usize]; 3],
            p_mismatch,
            em_probs: [1.0 - p_mismatch, p_mismatch],
            max_states,
            states,
            n_states: 0,
            fwd: [vec![0.0; ms], vec![0.0; ms], vec![0.0; ms]],
            bwd: [vec![0.0; ms], vec![0.0; ms], vec![0.0; ms]],
            fwd_sums: [0.0; 3],
            bwd_miss1: Vec::new(),
            bwd_miss2: Vec::new(),
            bwd_het1: Vec::new(),
            bwd_het2: Vec::new(),
            swap_haps: false,
            n_swaps: 0,
        }
    }

    fn ensure_capacity(&mut self, n_unph: i32, n_miss: i32) {
        let ms = self.max_states as usize;
        while (self.ref_alleles.len() as i32) < n_miss {
            self.ref_alleles.push(vec![0; ms]);
            self.bwd_miss1.push(vec![0.0; ms]);
            self.bwd_miss2.push(vec![0.0; ms]);
        }
        while (self.bwd_het1.len() as i32) < n_unph {
            self.bwd_het1.push(vec![0.0; ms]);
            self.bwd_het2.push(vec![0.0; ms]);
        }
    }

    fn initialize_bwd_fields(&mut self) {
        let ns = self.n_states as usize;
        let v = 1.0f32 / self.n_states as f32;
        self.bwd[0][..ns].fill(v);
        let (b0, rest) = self.bwd.split_at_mut(1);
        rest[0][..ns].copy_from_slice(&b0[0][..ns]);
        rest[1][..ns].copy_from_slice(&b0[0][..ns]);
    }

    fn bwd_alg(&mut self, mc: &MarkerCluster) {
        let sample_phase = mc.sample_phase();
        let mut miss_index = sample_phase.n_missing() + sample_phase.n_masked() - 1;
        let mut unph_index = sample_phase.n_unphased() - 1;
        self.initialize_bwd_fields();
        let ns = self.n_states as usize;
        let last_cluster = mc.n_clusters() - 1;
        if mc.is_missing_gt_or_masked_het(last_cluster) {
            let mi = miss_index as usize;
            self.bwd_miss1[mi][..ns].copy_from_slice(&self.bwd[0][..ns]);
            self.bwd_miss2[mi][..ns].copy_from_slice(&self.bwd[0][..ns]);
            miss_index -= 1;
        }
        for c in (0..last_cluster).rev() {
            self.bwd_step(mc, c);
            if mc.is_missing_gt_or_masked_het(c) {
                let mi = miss_index as usize;
                self.bwd_miss1[mi][..ns].copy_from_slice(&self.bwd[1][..ns]);
                self.bwd_miss2[mi][..ns].copy_from_slice(&self.bwd[2][..ns]);
                miss_index -= 1;
            }
            if mc.is_unphased_het(c + 1) {
                let ui = unph_index as usize;
                self.bwd_het1[ui][..ns].copy_from_slice(&self.bwd[1][..ns]);
                self.bwd_het2[ui][..ns].copy_from_slice(&self.bwd[2][..ns]);
                let (b0, rest) = self.bwd.split_at_mut(1);
                rest[0][..ns].copy_from_slice(&b0[0][..ns]); // bwd[1] = bwd[0]
                rest[1][..ns].copy_from_slice(&b0[0][..ns]); // bwd[2] = bwd[0]
                unph_index -= 1;
            }
        }
        debug_assert_eq!(miss_index, -1);
        debug_assert_eq!(unph_index, -1);
    }

    fn bwd_step(&mut self, mc: &MarkerCluster, cluster: i32) {
        let c_p1 = cluster + 1;
        let p_rec = mc.p_recomb().get(c_p1);
        let mut clust_em = (mc.cluster_end(c_p1) - mc.cluster_start(c_p1)) as f32 * self.p_mismatch;
        if clust_em >= 0.5 {
            clust_em = 0.5;
        }
        self.em_probs[1] = clust_em;
        self.em_probs[0] = 1.0 - clust_em;
        let ns = self.n_states;
        let cp = c_p1 as usize;
        HmmUpdater::bwd_update(
            &mut self.bwd[0],
            p_rec,
            &self.em_probs,
            &self.mismatches[0][cp],
            ns,
        );
        HmmUpdater::bwd_update(
            &mut self.bwd[1],
            p_rec,
            &self.em_probs,
            &self.mismatches[1][cp],
            ns,
        );
        HmmUpdater::bwd_update(
            &mut self.bwd[2],
            p_rec,
            &self.em_probs,
            &self.mismatches[2][cp],
            ns,
        );
    }

    fn initialize_fwd_fields(&mut self) {
        let ns = self.n_states as usize;
        let v = 1.0f32 / self.n_states as f32;
        self.fwd[0][..ns].fill(v);
        let (f0, rest) = self.fwd.split_at_mut(1);
        rest[0][..ns].copy_from_slice(&f0[0][..ns]);
        rest[1][..ns].copy_from_slice(&f0[0][..ns]);
        self.fwd_sums = [1.0, 1.0, 1.0];
    }

    fn fwd_alg(&mut self, mc: &mut MarkerCluster) {
        let mut miss_index = 0;
        let mut unph_het_index = 0;
        self.initialize_fwd_fields();
        let ns = self.n_states as usize;
        let n_clusters = mc.n_clusters();
        for c in 0..n_clusters {
            if mc.is_unphased_het(c) {
                self.phase_het(mc, unph_het_index, c);
                unph_het_index += 1;
                if self.swap_haps {
                    let swap_end = if unph_het_index < mc.unphased_het_clusters().size() {
                        mc.unphased_het_clusters().get(unph_het_index)
                    } else {
                        mc.n_clusters()
                    };
                    self.do_swap_haps(mc, c, swap_end);
                }
                let (f0, rest) = self.fwd.split_at_mut(1);
                rest[0][..ns].copy_from_slice(&f0[0][..ns]);
                rest[1][..ns].copy_from_slice(&f0[0][..ns]);
                self.fwd_sums[1] = self.fwd_sums[0];
                self.fwd_sums[2] = self.fwd_sums[0];
            }
            self.fwd_step(mc, c);
            if mc.is_missing_gt_or_masked_het(c) {
                self.impute_alleles(mc, c, miss_index);
                miss_index += 1;
            }
        }
    }

    fn fwd_step(&mut self, mc: &MarkerCluster, cluster: i32) {
        let p_rec = mc.p_recomb().get(cluster);
        let mut clust_em =
            (mc.cluster_end(cluster) - mc.cluster_start(cluster)) as f32 * self.p_mismatch;
        if clust_em >= 0.5 {
            clust_em = 0.5;
        }
        self.em_probs[1] = clust_em;
        self.em_probs[0] = 1.0 - clust_em;
        let ns = self.n_states;
        let cu = cluster as usize;
        self.fwd_sums[0] = HmmUpdater::fwd_update(
            &mut self.fwd[0],
            self.fwd_sums[0],
            p_rec,
            &self.em_probs,
            &self.mismatches[0][cu],
            ns,
        );
        self.fwd_sums[1] = HmmUpdater::fwd_update(
            &mut self.fwd[1],
            self.fwd_sums[1],
            p_rec,
            &self.em_probs,
            &self.mismatches[1][cu],
            ns,
        );
        self.fwd_sums[2] = HmmUpdater::fwd_update(
            &mut self.fwd[2],
            self.fwd_sums[2],
            p_rec,
            &self.em_probs,
            &self.mismatches[2][cu],
            ns,
        );
    }

    fn do_swap_haps(&mut self, mc: &mut MarkerCluster, start_clust: i32, end_clust: i32) {
        let (left, right) = self.mismatches.split_at_mut(2); // left=[0,1], right=[2]
        for c in start_clust..end_clust {
            std::mem::swap(&mut left[1][c as usize], &mut right[0][c as usize]);
        }
        let cs = mc.cluster_start(start_clust);
        let ce = mc.cluster_end(end_clust - 1);
        mc.sample_phase_mut().swap_haps(cs, ce);
    }

    fn phase_het(&mut self, mc: &mut MarkerCluster, unph_het_index: i32, cluster: i32) {
        let ns = self.n_states as usize;
        let idx = unph_het_index as usize;
        let mut p11 = 0.0f32;
        let mut p12 = 0.0f32;
        let mut p21 = 0.0f32;
        let mut p22 = 0.0f32;
        for k in 0..ns {
            let f1 = self.fwd[1][k];
            let f2 = self.fwd[2][k];
            let b1 = self.bwd_het1[idx][k];
            let b2 = self.bwd_het2[idx][k];
            p11 += f1 * b1;
            p12 += f1 * b2;
            p21 += f2 * b1;
            p22 += f2 * b2;
        }
        let num = p11 * p22;
        let den = p12 * p21;
        let last_swap_haps = self.swap_haps;
        self.swap_haps = num < den;
        if self.swap_haps != last_swap_haps {
            self.n_swaps += 1;
        }
        if !self.burnin
            && (num >= den * self.lr_threshold
                || (self.swap_haps && den >= num * self.lr_threshold))
        {
            mc.sample_phase_mut()
                .mark_unphased_het_cluster_as_phased(cluster);
        }
    }

    fn impute_alleles(&mut self, mc: &mut MarkerCluster, cluster: i32, miss_index: i32) {
        let ns = self.n_states as usize;
        let mi = miss_index as usize;
        let marker = mc.cluster_start(cluster);
        let n_alleles = self.markers.marker(marker).n_alleles() as usize;
        let mut al_freq1 = vec![0.0f32; n_alleles];
        let mut al_freq2 = vec![0.0f32; n_alleles];
        // stateProbs1/2 multiplied by fwd[1]/fwd[2]; if swapHaps, the two buffers swap roles.
        if self.swap_haps {
            for k in 0..ns {
                self.bwd_miss2[mi][k] *= self.fwd[1][k];
                self.bwd_miss1[mi][k] *= self.fwd[2][k];
            }
            for k in 0..ns {
                let a = self.ref_alleles[mi][k] as usize;
                al_freq1[a] += self.bwd_miss2[mi][k];
                al_freq2[a] += self.bwd_miss1[mi][k];
            }
        } else {
            for k in 0..ns {
                self.bwd_miss1[mi][k] *= self.fwd[1][k];
                self.bwd_miss2[mi][k] *= self.fwd[2][k];
            }
            for k in 0..ns {
                let a = self.ref_alleles[mi][k] as usize;
                al_freq1[a] += self.bwd_miss1[mi][k];
                al_freq2[a] += self.bwd_miss2[mi][k];
            }
        }
        let clust_type = mc.sample_phase().clust_type(cluster);
        if clust_type == ClustType::MissingGt {
            impute_missing_gt(mc.sample_phase_mut(), marker, &al_freq1, &al_freq2);
        } else if clust_type == ClustType::MaskedHet {
            self.impute_masked_het(mc, cluster, marker, &al_freq1, &al_freq2);
        }
    }

    fn impute_masked_het(
        &self,
        mc: &mut MarkerCluster,
        cluster: i32,
        marker: i32,
        al_freq1: &[f32],
        al_freq2: &[f32],
    ) {
        let lr_threshold = self.lr_threshold;
        let sample_phase = mc.sample_phase_mut();
        let a1 = sample_phase.allele1(marker);
        let a2 = sample_phase.allele2(marker);
        debug_assert_ne!(a1, a2);
        let p_no_switch = al_freq1[a1 as usize] * al_freq2[a2 as usize];
        let p_switch = al_freq1[a2 as usize] * al_freq2[a1 as usize];
        if p_switch > p_no_switch {
            sample_phase.set_allele1(marker, a2);
            sample_phase.set_allele2(marker, a1);
            if p_switch >= lr_threshold * p_no_switch {
                sample_phase.mark_masked_het_cluster_as_phased(cluster);
            }
        } else if p_no_switch >= lr_threshold * p_switch {
            sample_phase.mark_masked_het_cluster_as_phased(cluster);
        }
    }
}

fn impute_missing_gt(
    sample_phase: &mut SamplePhase,
    marker: i32,
    al_freq1: &[f32],
    al_freq2: &[f32],
) {
    let mut a1 = 0usize;
    let mut a2 = 0usize;
    for j in 1..al_freq1.len() {
        if al_freq1[j] > al_freq1[a1] {
            a1 = j;
        }
        if al_freq2[j] > al_freq2[a2] {
            a2 = j;
        }
    }
    sample_phase.set_allele1(marker, a1 as i32);
    sample_phase.set_allele2(marker, a2 as i32);
}

impl PhaseBaum for PhaseBaum2<'_> {
    fn n_targ_samples(&self) -> i32 {
        self.phase_ibs.phase_data().fpd().targ_gt().n_samples()
    }

    fn phase(&mut self, sample: i32) {
        let phase_ibs = self.phase_ibs;
        let phase_data = phase_ibs.phase_data();
        let est_phase = phase_data.est_phase();
        let mut sample_phase = est_phase.get(sample);
        if self.mask_trailing_hets {
            sample_phase.mask_trailing_unphased_hets();
        }
        let n_unph_hets = sample_phase.n_unphased();
        let n_masked_hets = sample_phase.n_masked();
        let n_missing_or_masked = sample_phase.n_missing() + n_masked_hets;
        if n_missing_or_masked > 0 || n_unph_hets > 0 {
            self.n_swaps = 0;
            self.swap_haps = false;
            // Store the masked clone first so the MarkerCluster sees it (Java aliases the
            // shared SamplePhase object).
            est_phase.set(sample, sample_phase);
            let mut mc = MarkerCluster::new(phase_data, sample);
            self.ensure_capacity(n_unph_hets, n_missing_or_masked);
            self.n_states =
                self.states
                    .ibs_states(&mc, &mut self.ref_alleles, &mut self.mismatches);
            self.bwd_alg(&mc);
            self.fwd_alg(&mut mc);
            est_phase.set(sample, mc.into_sample_phase());
            SwapRate::increment(n_unph_hets, self.n_swaps);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::main_pkg::{Par, Pedigree};
    use crate::vcf::{
        BasicGT, BasicGTRec, GTRec, GeneticMap, MarkerIndices, MarkerParser, PositionMap,
        VcfHeader, VcfRecGTParser, Window, GT, HEADER_PREFIX,
    };
    use std::io::Write;
    use std::rc::Rc;

    use super::super::{FixedPhaseData, PhaseData};

    fn fpd() -> FixedPhaseData {
        let gt = std::env::temp_dir().join("beagle_rs_pb2_gt.vcf");
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
        // include a het sample and a missing genotype to exercise the HMM branches
        let gts = [
            "0|0\t0|1\t1|1\t0|1\t0|0\t1|0",
            "0|1\t0|0\t1|0\t1|1\t0|1\t0|0",
            "1|1\t0|1\t.|.\t0|1\t1|0\t0|1",
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
    fn phase_runs_and_preserves_genotypes() {
        let pd = PhaseData::new(Rc::new(fpd()), 99999);
        let n_markers = pd.fpd().stage1_targ_gt().n_markers();
        let n_samples = pd.fpd().stage1_targ_gt().n_samples();
        // snapshot the unordered genotype (allele set) per sample/marker before phasing
        let before: Vec<Vec<(i32, i32)>> = (0..n_samples)
            .map(|s| {
                (0..n_markers)
                    .map(|m| {
                        let a = pd.est_phase().get(s).allele1(m);
                        let b = pd.est_phase().get(s).allele2(m);
                        if a <= b {
                            (a, b)
                        } else {
                            (b, a)
                        }
                    })
                    .collect()
            })
            .collect();

        let phase_ibs = PbwtPhaseIbs::new(&pd, pd.coded_steps(), false);
        let mut baum = PhaseBaum2::new(&phase_ibs);
        assert_eq!(baum.n_targ_samples(), n_samples);
        for s in 0..baum.n_targ_samples() {
            baum.phase(s);
        }

        // after phasing, each non-missing genotype's allele set is unchanged (phase may flip);
        // the previously-missing genotype (sample 2, marker 2) is now imputed (non-negative).
        for s in 0..n_samples {
            for m in 0..n_markers {
                let a = pd.est_phase().get(s).allele1(m);
                let b = pd.est_phase().get(s).allele2(m);
                assert!(a >= 0 && b >= 0, "s={s} m={m} not imputed");
                let set = if a <= b { (a, b) } else { (b, a) };
                let was = before[s as usize][m as usize];
                if was.0 >= 0 && was.1 >= 0 {
                    assert_eq!(set, was, "s={s} m={m} genotype changed");
                }
            }
        }
    }
}
