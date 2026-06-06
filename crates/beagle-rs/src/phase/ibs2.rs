//! Port of `phase/Ibs2.java` — the IBS2 segments each target sample shares with another
//! target sample (merge gaps ≤ 4 cM, extend by phase-consistency, keep segments ≥ 2 cM).

use crate::blbutil::{DoubleArray, FloatArray};
use crate::vcf::{MarkerMap, GT};

use super::{Ibs2Markers, Ibs2Sets, SampleSeg};

const MIN_IBS2_CM: f32 = 2.0;
const MAX_IBD_GAP_CM: f32 = 4.0;

/// Port of `phase/Ibs2.java`.
pub struct Ibs2 {
    n_markers: i32,
    sample_segs: Vec<Vec<SampleSeg>>, // [targ sample][segment]
}

fn are_phase_consistent(a1: i32, a2: i32, b1: i32, b2: i32) -> bool {
    (a1 < 0 || b1 < 0 || a1 == b1) && (a2 < 0 || b2 < 0 || a2 == b2)
}

fn ibs2(targ_gt: &dyn GT, m: i32, s1: i32, s2: i32) -> bool {
    let hap1 = s1 << 1;
    let hap2 = s2 << 1;
    let a1 = targ_gt.allele(m, hap1);
    let a2 = targ_gt.allele(m, hap1 | 0b1);
    let b1 = targ_gt.allele(m, hap2);
    let b2 = targ_gt.allele(m, hap2 | 0b1);
    are_phase_consistent(a1, a2, b1, b2) || are_phase_consistent(a1, a2, b2, b1)
}

fn gap_cm(prev: &SampleSeg, next: &SampleSeg, gen_pos: &DoubleArray) -> f64 {
    use crate::beagleutil::IntInterval;
    gen_pos.get(next.start()) - gen_pos.get(prev.incl_end())
}

fn merge_segments(list: Vec<SampleSeg>, gen_pos: &DoubleArray) -> Vec<SampleSeg> {
    use crate::beagleutil::IntInterval;
    if list.len() < 2 {
        return list;
    }
    let mut merged = Vec::new();
    let mut prev = list[0];
    for next in &list[1..] {
        if prev.sample() == next.sample() && gap_cm(&prev, next, gen_pos) <= MAX_IBD_GAP_CM as f64 {
            debug_assert!(prev.incl_end() <= next.incl_end());
            prev = SampleSeg::new(prev.sample(), prev.start(), next.incl_end());
        } else {
            merged.push(prev);
            prev = *next;
        }
    }
    merged.push(prev);
    merged
}

fn extend(targ_gt: &dyn GT, sample: i32, ss: &SampleSeg) -> SampleSeg {
    use crate::beagleutil::IntInterval;
    let n_markers = targ_gt.n_markers();
    let sample2 = ss.sample();
    let mut incl_start = ss.start();
    let mut excl_end = ss.incl_end() + 1;
    while incl_start > 0 && ibs2(targ_gt, incl_start - 1, sample, sample2) {
        incl_start -= 1;
    }
    while excl_end < n_markers && ibs2(targ_gt, excl_end, sample, sample2) {
        excl_end += 1;
    }
    SampleSeg::new(sample2, incl_start, excl_end - 1)
}

fn apply_length_filter(list: Vec<SampleSeg>, gen_pos: &DoubleArray) -> Vec<SampleSeg> {
    use crate::beagleutil::IntInterval;
    list.into_iter()
        .filter(|ss| (gen_pos.get(ss.incl_end()) - gen_pos.get(ss.start())) >= MIN_IBS2_CM as f64)
        .collect()
}

fn ibs2_segments(
    targ_gt: &dyn GT,
    gen_pos: &DoubleArray,
    ibs2_sets: &Ibs2Sets,
    sample: i32,
) -> Vec<SampleSeg> {
    let mut seg_list = ibs2_sets.seg_list(sample);
    seg_list.sort_by(SampleSeg::sample_cmp);
    let seg_list = merge_segments(seg_list, gen_pos);
    let seg_list: Vec<SampleSeg> = seg_list
        .iter()
        .map(|ss| extend(targ_gt, sample, ss))
        .collect();
    let seg_list = merge_segments(seg_list, gen_pos);
    apply_length_filter(seg_list, gen_pos)
}

impl Ibs2 {
    /// `new Ibs2(GT targGT, MarkerMap map, FloatArray maf)`.
    pub fn new(targ_gt: &dyn GT, map: &MarkerMap, maf: &FloatArray) -> Self {
        let ibs2_markers = Ibs2Markers::new(targ_gt, map, maf);
        let ibs2_sets = Ibs2Sets::new(targ_gt, &ibs2_markers);
        let gen_pos = map.gen_pos();
        let sample_segs: Vec<Vec<SampleSeg>> = (0..targ_gt.n_samples())
            .map(|s| ibs2_segments(targ_gt, gen_pos, &ibs2_sets, s))
            .collect();
        Ibs2 {
            n_markers: targ_gt.n_markers(),
            sample_segs,
        }
    }

    /// `nMarkers()`.
    pub fn n_markers(&self) -> i32 {
        self.n_markers
    }

    /// `nTargSamples()`.
    pub fn n_targ_samples(&self) -> i32 {
        self.sample_segs.len() as i32
    }

    /// `nIbs2Segments(int targSample)`.
    pub fn n_ibs2_segments(&self, targ_sample: i32) -> i32 {
        self.sample_segs[targ_sample as usize].len() as i32
    }

    /// `areIbs2(int targSample, int otherSample, int marker)`.
    pub fn are_ibs2_at(&self, targ_sample: i32, other_sample: i32, marker: i32) -> bool {
        use crate::beagleutil::IntInterval;
        assert!(
            (targ_sample as usize) < self.sample_segs.len(),
            "{targ_sample}"
        );
        assert!(other_sample >= 0, "{other_sample}");
        assert!(marker >= 0 && marker < self.n_markers, "{marker}");
        if targ_sample == other_sample {
            return true;
        }
        if !self.sample_segs[targ_sample as usize].is_empty()
            && (other_sample as usize) < self.sample_segs.len()
        {
            for ss in &self.sample_segs[targ_sample as usize] {
                if ss.sample() == other_sample && ss.start() <= marker && marker <= ss.incl_end() {
                    return true;
                }
            }
        }
        false
    }

    /// `areIbs2(int targSample, int otherSample, int start, int inclEnd)`.
    pub fn are_ibs2_in(
        &self,
        targ_sample: i32,
        other_sample: i32,
        start: i32,
        incl_end: i32,
    ) -> bool {
        use crate::beagleutil::IntInterval;
        assert!(start <= incl_end, "{start}");
        let same_sample = targ_sample == other_sample;
        if self.sample_segs[targ_sample as usize].is_empty() || same_sample {
            return same_sample;
        }
        for ss in &self.sample_segs[targ_sample as usize] {
            if ss.sample() == other_sample && start <= ss.incl_end() && ss.start() <= incl_end {
                return true;
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::{
        BasicGT, BasicGTRec, GTRec, MarkerParser, PositionMap, VcfHeader, VcfRecGTParser,
        HEADER_PREFIX,
    };
    use std::rc::Rc;

    #[test]
    fn ibs2_between_identical_het_samples_spans_window() {
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..3 {
            hdr.push_str(&format!("\tS{s}"));
        }
        let header = VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &[true; 3],
        );
        let mp = MarkerParser::new(true, true, true, true);
        let n = 60;
        let recs: Vec<Rc<dyn GTRec>> = (0..n)
            .map(|i| {
                let pos = 100_000 + i * 50_000;
                let line = format!("chr1\t{pos}\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t0|1\t0|0");
                let p = VcfRecGTParser::new(&header, &line, &mp);
                Rc::new(BasicGTRec::from_parser(&p)) as Rc<dyn GTRec>
            })
            .collect();
        let gt = BasicGT::new(recs);
        let map = MarkerMap::create_min_dist(&PositionMap::new(1e-6), 0.0, gt.markers());
        let maf = FloatArray::from_doubles(&vec![1.0 / 3.0; n as usize]);
        let ibs2 = Ibs2::new(&gt, &map, &maf);

        assert_eq!(ibs2.n_targ_samples(), 3);
        assert_eq!(ibs2.n_markers(), n);
        // a sample is trivially IBS2 with itself
        assert!(ibs2.are_ibs2_at(0, 0, 30));
        // samples 0 & 1 (identical het) are IBS2 across the segment
        assert!(ibs2.n_ibs2_segments(0) >= 1);
        assert!(ibs2.are_ibs2_at(0, 1, 30));
        assert!(ibs2.are_ibs2_in(0, 1, 10, 40));
        // sample 2 (homozygous, different) is not IBS2 with sample 0
        assert!(!ibs2.are_ibs2_at(0, 2, 30));
    }
}
