//! Port of `phase/FixedPhaseData.java` — immutable per-window phasing input: the marker
//! map, target/reference genotypes, the stage-1 (high-frequency) marker subset, MAF
//! estimates, carrier lists, and IBS2 data.

use std::rc::Rc;

use crate::blbutil::{DoubleArray, FloatArray, Utilities};
use crate::ints::{java_binary_search, IntArray, WrappedIntArray};
use crate::jdk::Random;
use crate::main_pkg::{Par, Pedigree};
use crate::vcf::{
    mean_single_base_gen_dist, CarrierList, GeneticMap, MarkerMap, Markers, RefGT, SplicedGT,
    Steps, Window, XRefGT, GT,
};

use super::Ibs2;

const MAX_HIFREQ_PROP: f32 = 0.75;

/// Port of `phase/FixedPhaseData.java`.
pub struct FixedPhaseData {
    par: Par,
    ped: Pedigree,
    window: i32,
    map: Rc<MarkerMap>,
    stage1_steps: Steps,
    targ_gt: Rc<dyn GT>,
    restricted_ref_gt: Option<RefGT>,
    overlap: i32,
    stage1_map: Rc<MarkerMap>,
    ibs_step: f32,
    stage1_targ_gt: Rc<dyn GT>,
    stage1_ref_gt: Option<RefGT>,
    stage1_xref_gt: Option<XRefGT>,
    stage1_maf: FloatArray,
    stage1_overlap: i32,
    stage1_ibs2: Ibs2,
    n_haps: i32,
    carriers: Vec<Vec<CarrierList>>,
    stage1_to2: WrappedIntArray,
    prev_stage1_marker: Vec<i32>,
    prev_stage1_wt: Vec<f32>,
}

fn marker_map(gen_map: &dyn GeneticMap, markers: &Markers) -> MarkerMap {
    let mean_gen_diff = mean_single_base_gen_dist(gen_map, markers);
    MarkerMap::create_min_dist(gen_map, mean_gen_diff, markers)
}

fn median_diff(da: &DoubleArray) -> f32 {
    let mut diffs: Vec<f64> = (1..da.size()).map(|j| da.get(j) - da.get(j - 1)).collect();
    diffs.sort_by(|a, b| a.total_cmp(b));
    let n = diffs.len();
    0.5f32 * ((diffs[(n - 1) >> 1] + diffs[n >> 1]) as f32)
}

fn n_haps(window: &Window) -> i32 {
    let n_ref_haps = window.ref_gt().map_or(0, |r| r.n_haps());
    window.targ_gt().n_haps() + n_ref_haps
}

fn carriers(par: &Par, window: &Window) -> Vec<Vec<CarrierList>> {
    let n_ref_samples = window.ref_gt().map_or(0, |r| r.n_samples());
    let n_samples = window.targ_gt().n_samples() + n_ref_samples;
    let max_carriers = 3.max(((n_samples as f32 * par.rare()) as f64).floor() as i32);
    window.carriers(max_carriers)
}

fn hi_freq_indices(carriers: &[Vec<CarrierList>]) -> Vec<i32> {
    (0..carriers.len() as i32)
        .filter(|&m| {
            carriers[m as usize]
                .iter()
                .filter(|c| matches!(c, CarrierList::HighFreq))
                .count()
                > 1
        })
        .collect()
}

fn ignore_low_freq_carriers(carriers: &mut [Vec<CarrierList>]) {
    for row in carriers.iter_mut() {
        for c in row.iter_mut() {
            *c = CarrierList::HighFreq;
        }
    }
}

fn rand_haps(gt: &dyn GT, max_haps: i32, rand: &mut Random) -> Vec<i32> {
    let n_haps = gt.n_haps();
    let mut ia: Vec<i32> = (0..n_haps).collect();
    if n_haps > max_haps {
        Utilities::shuffle_n(&mut ia, max_haps, rand);
        ia.truncate(max_haps as usize);
        ia.sort_unstable();
    }
    ia
}

fn maf_at(
    gt: &dyn GT,
    opt_ref_gt: Option<&RefGT>,
    targ_haps: &[i32],
    ref_haps: &[i32],
    m: i32,
) -> f64 {
    let mut mod_cnts = vec![0i32; (gt.marker(m).n_alleles() + 1) as usize];
    for &h in targ_haps {
        mod_cnts[(gt.allele(m, h) + 1) as usize] += 1;
    }
    if let Some(ref_gt) = opt_ref_gt {
        if !ref_haps.is_empty() {
            for &h in ref_haps {
                mod_cnts[(ref_gt.allele(m, h) + 1) as usize] += 1;
            }
        }
    }
    mod_cnts[0] = 0; // zero out missing count
    mod_cnts.sort_unstable();
    let den: i32 = mod_cnts[1..].iter().sum();
    if den == 0 {
        0.0
    } else {
        mod_cnts[mod_cnts.len() - 2] as f64 / den as f64
    }
}

fn maf(opt_ref_gt: Option<&RefGT>, targ_gt: &dyn GT, max_haps: i32, seed: i64) -> FloatArray {
    let mut rand = Random::new(seed);
    let targ_haps = rand_haps(targ_gt, max_haps, &mut rand);
    let ref_haps = if (targ_haps.len() as i32) < max_haps {
        match opt_ref_gt {
            Some(r) => rand_haps(r, max_haps - targ_haps.len() as i32, &mut rand),
            None => Vec::new(),
        }
    } else {
        Vec::new()
    };
    let maf: Vec<f64> = (0..targ_gt.n_markers())
        .map(|m| maf_at(targ_gt, opt_ref_gt, &targ_haps, &ref_haps, m))
        .collect();
    FloatArray::from_doubles(&maf)
}

fn stage1_targ_overlap(phased_overlap: Option<&Rc<dyn GT>>, hi_freq_mkrs: &[i32]) -> i32 {
    match phased_overlap {
        None => 0,
        Some(po) => {
            let ins_pt =
                java_binary_search(hi_freq_mkrs, 0, hi_freq_mkrs.len() as i32, po.n_markers());
            if ins_pt < 0 {
                -ins_pt - 1
            } else {
                ins_pt
            }
        }
    }
}

fn prev_stage1_marker(n_markers: i32, stage1_indices: &WrappedIntArray) -> Vec<i32> {
    let mut mkr_a = vec![0i32; n_markers as usize];
    let n_hi_freq = stage1_indices.size();
    let mut start = stage1_indices.get(1);
    for j in 2..n_hi_freq {
        let end = stage1_indices.get(j);
        for m in start..end {
            mkr_a[m as usize] = j - 1;
        }
        start = end;
    }
    for m in start..n_markers {
        mkr_a[m as usize] = n_hi_freq - 1;
    }
    mkr_a
}

fn prev_wt(map: &MarkerMap, marker_indices: &WrappedIntArray) -> Vec<f32> {
    let gen_pos = map.gen_pos();
    let mut prev_wt = vec![0.0f32; gen_pos.size() as usize];
    for w in prev_wt.iter_mut().take(marker_indices.get(0) as usize) {
        *w = 1.0;
    }
    let mut start = marker_indices.get(0);
    for j in 1..marker_indices.size() {
        let end = marker_indices.get(j);
        let pos_a = gen_pos.get(start);
        let pos_b = gen_pos.get(end);
        let d = pos_b - pos_a;
        prev_wt[start as usize] = 1.0;
        for m in (start + 1)..end {
            prev_wt[m as usize] = ((pos_b - gen_pos.get(m)) / d) as f32;
        }
        start = end;
    }
    for m in start..gen_pos.size() {
        prev_wt[m as usize] = 1.0;
    }
    prev_wt
}

fn check_data(window: &Window, phased_overlap: Option<&Rc<dyn GT>>) {
    if let Some(po) = phased_overlap {
        let targ = window.targ_gt();
        assert!(po.is_phased(), "unphased");
        assert!(po.samples() == targ.samples(), "inconsistent data");
        assert!(po.n_markers() <= targ.n_markers(), "inconsistent data");
        for j in 0..po.n_markers() {
            assert!(po.marker(j) == targ.marker(j), "inconsistent data");
        }
    }
}

impl FixedPhaseData {
    /// `new FixedPhaseData(Par par, Pedigree ped, Window window, GT phasedOverlap)`.
    pub fn new(
        par: &Par,
        ped: &Pedigree,
        window: &Window,
        phased_overlap: Option<Rc<dyn GT>>,
    ) -> Self {
        check_data(window, phased_overlap.as_ref());
        let n_targ_markers = window.targ_gt().n_markers();
        let window_idx = window.window_index();

        let map = Rc::new(marker_map(window.gen_map(), window.targ_gt().markers()));
        let targ_base: Rc<dyn GT> = Rc::new(window.targ_gt().clone());
        let targ_gt: Rc<dyn GT> = match &phased_overlap {
            None => targ_base,
            Some(po) => Rc::new(SplicedGT::new(po.clone(), targ_base)),
        };
        let restricted_ref_gt = window.restrict_ref_gt().cloned();
        let overlap = phased_overlap.as_ref().map_or(0, |po| po.n_markers());
        let n_haps_val = n_haps(window);

        let mut rare_carriers = carriers(par, window);
        let mut hi_freq_ind = hi_freq_indices(&rare_carriers);

        let (
            carriers_field,
            stage1_map,
            ibs_step,
            stage1_steps,
            stage1_targ_gt,
            stage1_ref_gt,
            stage1_overlap,
            stage1_to2,
            prev_stage1_marker_v,
            prev_stage1_wt_v,
        );

        if (hi_freq_ind.len() as i32) < 2
            || hi_freq_ind.len() as f32 > MAX_HIFREQ_PROP * n_targ_markers as f32
        {
            hi_freq_ind = (0..n_targ_markers).collect();
            ignore_low_freq_carriers(&mut rare_carriers);
            carriers_field = rare_carriers;
            let s1_map = map.clone();
            ibs_step = par.step_scale() * median_diff(s1_map.gen_pos());
            stage1_steps = Steps::new(s1_map.clone(), ibs_step);
            stage1_map = s1_map;
            stage1_targ_gt = targ_gt.clone();
            stage1_ref_gt = restricted_ref_gt.clone();
            stage1_overlap = overlap;
            stage1_to2 = WrappedIntArray::from_slice(&hi_freq_ind);
            prev_stage1_marker_v = (0..targ_gt.n_markers()).collect();
            prev_stage1_wt_v = vec![1.0f32; targ_gt.n_markers() as usize];
        } else {
            let hi_freq_markers = targ_gt.markers().restrict_indices(&hi_freq_ind);
            carriers_field = rare_carriers;
            let s1_map = Rc::new(map.restrict(&hi_freq_ind));
            ibs_step = par.step_scale() * median_diff(s1_map.gen_pos());
            stage1_steps = Steps::new(s1_map.clone(), ibs_step);
            stage1_targ_gt = targ_gt.clone().restrict(&hi_freq_markers, &hi_freq_ind);
            stage1_ref_gt = restricted_ref_gt
                .as_ref()
                .map(|r| r.restrict_to_ref(&hi_freq_markers, &hi_freq_ind));
            stage1_overlap = stage1_targ_overlap(phased_overlap.as_ref(), &hi_freq_ind);
            stage1_to2 = WrappedIntArray::from_slice(&hi_freq_ind);
            prev_stage1_marker_v = prev_stage1_marker(targ_gt.n_markers(), &stage1_to2);
            prev_stage1_wt_v = prev_wt(&map, &stage1_to2);
            stage1_map = s1_map;
        }

        let stage1_xref_gt = stage1_ref_gt
            .as_ref()
            .map(|r| XRefGT::from_phased_gt(r, par.nthreads()));

        let max_maf_haps = 10000;
        let stage1_maf = maf(
            stage1_ref_gt.as_ref(),
            stage1_targ_gt.as_ref(),
            max_maf_haps,
            par.seed(),
        );
        let stage1_ibs2 = Ibs2::new(stage1_targ_gt.as_ref(), &stage1_map, &stage1_maf);

        FixedPhaseData {
            par: par.clone(),
            ped: ped.clone(),
            window: window_idx,
            map,
            stage1_steps,
            targ_gt,
            restricted_ref_gt,
            overlap,
            stage1_map,
            ibs_step,
            stage1_targ_gt,
            stage1_ref_gt,
            stage1_xref_gt,
            stage1_maf,
            stage1_overlap,
            stage1_ibs2,
            n_haps: n_haps_val,
            carriers: carriers_field,
            stage1_to2,
            prev_stage1_marker: prev_stage1_marker_v,
            prev_stage1_wt: prev_stage1_wt_v,
        }
    }

    /// `par()`.
    pub fn par(&self) -> &Par {
        &self.par
    }
    /// `window()`.
    pub fn window(&self) -> i32 {
        self.window
    }
    /// `ped()`.
    pub fn ped(&self) -> &Pedigree {
        &self.ped
    }
    /// `map()`.
    pub fn map(&self) -> &MarkerMap {
        &self.map
    }
    /// `restrictedRefGT()`.
    pub fn restricted_ref_gt(&self) -> Option<&RefGT> {
        self.restricted_ref_gt.as_ref()
    }
    /// `targGT()`.
    pub fn targ_gt(&self) -> &Rc<dyn GT> {
        &self.targ_gt
    }
    /// `overlap()`.
    pub fn overlap(&self) -> i32 {
        self.overlap
    }
    /// `stage1Map()`.
    pub fn stage1_map(&self) -> &MarkerMap {
        &self.stage1_map
    }
    /// `ibsStep()`.
    pub fn ibs_step(&self) -> f32 {
        self.ibs_step
    }
    /// `stage1Steps()`.
    pub fn stage1_steps(&self) -> &Steps {
        &self.stage1_steps
    }
    /// `stage1RefGT()`.
    pub fn stage1_ref_gt(&self) -> Option<&RefGT> {
        self.stage1_ref_gt.as_ref()
    }
    /// `stage1XRefGT()`.
    pub fn stage1_xref_gt(&self) -> Option<&XRefGT> {
        self.stage1_xref_gt.as_ref()
    }
    /// `stage1TargGT()`.
    pub fn stage1_targ_gt(&self) -> &Rc<dyn GT> {
        &self.stage1_targ_gt
    }
    /// `stage1Maf()`.
    pub fn stage1_maf(&self) -> &FloatArray {
        &self.stage1_maf
    }
    /// `stage1Overlap()`.
    pub fn stage1_overlap(&self) -> i32 {
        self.stage1_overlap
    }
    /// `nHaps()`.
    pub fn n_haps(&self) -> i32 {
        self.n_haps
    }
    /// `stage1To2()`.
    pub fn stage1_to2(&self) -> &WrappedIntArray {
        &self.stage1_to2
    }
    /// `stage1Ibs2()`.
    pub fn stage1_ibs2(&self) -> &Ibs2 {
        &self.stage1_ibs2
    }
    /// `carriers(int marker, int allele)`.
    pub fn carriers(&self, marker: i32, allele: i32) -> &CarrierList {
        &self.carriers[marker as usize][allele as usize]
    }
    /// `isLowFreq(int marker, int allele)`.
    pub fn is_low_freq(&self, marker: i32, allele: i32) -> bool {
        !matches!(
            self.carriers[marker as usize][allele as usize],
            CarrierList::HighFreq
        )
    }
    /// `prevStage1Marker(int marker)`.
    pub fn prev_stage1_marker(&self, marker: i32) -> i32 {
        self.prev_stage1_marker[marker as usize]
    }
    /// `prevStage1Wt(int marker)`.
    pub fn prev_stage1_wt(&self, marker: i32) -> f32 {
        self.prev_stage1_wt[marker as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::{
        BasicGT, BasicGTRec, GTRec, MarkerIndices, MarkerParser, PositionMap, VcfHeader,
        VcfRecGTParser, HEADER_PREFIX,
    };
    use std::io::Write;

    fn make_par() -> Par {
        let gt = std::env::temp_dir().join("beagle_rs_fpd_gt.vcf");
        std::fs::File::create(&gt).unwrap().write_all(b"x").unwrap();
        Par::new(&[
            format!("gt={}", gt.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
            "seed=99999".to_string(),
        ])
    }

    fn window(lines: &[&str], n_dip: usize) -> Window {
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..n_dip {
            hdr.push_str(&format!("\tS{s}"));
        }
        let h = VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &vec![true; n_dip],
        );
        let mp = MarkerParser::new(true, true, true, true);
        let recs: Vec<Rc<dyn GTRec>> = lines
            .iter()
            .map(|l| {
                Rc::new(BasicGTRec::from_parser(&VcfRecGTParser::new(&h, l, &mp))) as Rc<dyn GTRec>
            })
            .collect();
        let targ = BasicGT::new(recs);
        let n = targ.n_markers();
        let indices = MarkerIndices::from_counts(0, n, n);
        let gen_map: Rc<dyn GeneticMap> = Rc::new(PositionMap::new(1e-6));
        Window::new(gen_map, 1, true, indices, None, targ)
    }

    #[test]
    fn builds_fixed_phase_data_target_only() {
        let par = make_par();
        let w = window(
            &[
                "chr1\t1000000\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0",
                "chr1\t2000000\t.\tG\tT\t.\tPASS\t.\tGT\t0|0\t0|1",
                "chr1\t3000000\t.\tA\tG\t.\tPASS\t.\tGT\t1|1\t0|0",
                "chr1\t4000000\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\t1|1",
            ],
            2,
        );
        let ped = Pedigree::new(w.targ_gt().samples().clone(), None);
        let fpd = FixedPhaseData::new(&par, &ped, &w, None);
        assert_eq!(fpd.window(), 1);
        assert_eq!(fpd.targ_gt().n_markers(), 4);
        assert_eq!(fpd.n_haps(), 4); // 2 diploid samples, no ref
        assert!(fpd.restricted_ref_gt().is_none());
        assert!(fpd.stage1_ref_gt().is_none());
        // small window -> branch A: stage1 == full set
        assert_eq!(fpd.stage1_targ_gt().n_markers(), 4);
        assert_eq!(fpd.stage1_to2().size(), 4);
        assert_eq!(fpd.stage1_maf().size(), 4);
        assert_eq!(fpd.stage1_ibs2().n_markers(), 4);
        assert_eq!(fpd.prev_stage1_marker(3), 3);
        assert_eq!(fpd.prev_stage1_wt(0), 1.0);
        assert_eq!(fpd.par().seed(), 99999);
    }
}
