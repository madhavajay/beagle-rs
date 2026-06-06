//! Port of `phase/Ibs2Markers.java` — selects the markers and step intervals used to
//! detect IBS2 segments (common, low-missingness markers spaced ≥ a minimum cM apart).

use crate::blbutil::FloatArray;
use crate::ints::{IntList, WrappedIntArray};
use crate::vcf::{MarkerMap, GT};

const MAX_MISS_FREQ: f32 = 0.1;
const MIN_MINOR_FREQ: f32 = 0.1;
const MIN_MARKER_CNT: i32 = 50;
const MIN_INTERMARKER_CM: f64 = 0.02;

/// Port of `phase/Ibs2Markers.java`.
pub struct Ibs2Markers {
    use_marker: Vec<bool>,
    step_starts: WrappedIntArray,
}

fn use_marker(targ_gt: &dyn GT, m: i32, maf: &FloatArray, max_miss_cnt: i32) -> bool {
    if maf.get(m) >= MIN_MINOR_FREQ {
        let mut miss_cnt = 0;
        for h in 0..targ_gt.n_haps() {
            if targ_gt.allele(m, h) < 0 {
                miss_cnt += 1;
            }
        }
        miss_cnt <= max_miss_cnt
    } else {
        false
    }
}

fn next_start(gen_pos: &crate::blbutil::DoubleArray, start: i32, use_marker0: &mut [bool]) -> i32 {
    let mut cm_pos = gen_pos.get(start);
    let mut min_cm_pos = cm_pos + MIN_INTERMARKER_CM;
    let mut next_start = start + 1;
    let mut mkr_cnt = 0;
    while (next_start as usize) < use_marker0.len() && mkr_cnt < MIN_MARKER_CNT {
        if use_marker0[next_start as usize] {
            cm_pos = gen_pos.get(next_start);
            if cm_pos < min_cm_pos {
                use_marker0[next_start as usize] = false;
            } else {
                mkr_cnt += 1;
                min_cm_pos = cm_pos + MIN_INTERMARKER_CM;
            }
        }
        next_start += 1;
    }
    next_start
}

fn step_starts(use_marker0: &mut [bool], map: &MarkerMap) -> WrappedIntArray {
    let gen_pos = map.gen_pos();
    let n_markers = gen_pos.size();
    let mut indices = IntList::with_capacity(gen_pos.size() >> 6);
    let mut last_start = 0;
    let mut next = next_start(gen_pos, last_start, use_marker0);
    // combines the last two steps (the final start is not added)
    while next < n_markers {
        indices.add(last_start);
        last_start = next;
        next = next_start(gen_pos, next, use_marker0);
    }
    WrappedIntArray::from_list(&indices)
}

impl Ibs2Markers {
    /// `new Ibs2Markers(GT targGT, MarkerMap map, FloatArray maf)`.
    pub fn new(targ_gt: &dyn GT, map: &MarkerMap, maf: &FloatArray) -> Self {
        assert!(
            map.gen_pos().size() == targ_gt.n_markers(),
            "{}",
            map.gen_pos().size()
        );
        assert!(maf.size() == targ_gt.n_markers(), "{}", maf.size());
        let max_miss = (MAX_MISS_FREQ * targ_gt.n_haps() as f32).ceil() as i32;
        let mut use_marker0: Vec<bool> = (0..targ_gt.n_markers())
            .map(|m| use_marker(targ_gt, m, maf, max_miss))
            .collect();
        let step_starts = step_starts(&mut use_marker0, map);
        Ibs2Markers {
            use_marker: use_marker0,
            step_starts,
        }
    }

    /// `nMarkers()`.
    pub fn n_markers(&self) -> i32 {
        self.use_marker.len() as i32
    }

    /// `markers(int start, int end)` — the used marker indices in `[start, end)`.
    pub fn markers(&self, start: i32, end: i32) -> Vec<i32> {
        assert!(start >= 0 && end <= self.use_marker.len() as i32, "{start}");
        assert!(end <= self.use_marker.len() as i32, "{end}");
        let mut markers = Vec::with_capacity((end - start) as usize);
        for m in start..end {
            if self.use_marker[m as usize] {
                markers.push(m);
            }
        }
        markers
    }

    /// `stepStarts()` — the first marker index of each step.
    pub fn step_starts(&self) -> &WrappedIntArray {
        &self.step_starts
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ints::IntArray;
    use crate::vcf::{
        BasicGT, BasicGTRec, GTRec, MarkerParser, PositionMap, VcfHeader, VcfRecGTParser,
        HEADER_PREFIX,
    };
    use std::rc::Rc;

    fn header(n_dip: usize) -> VcfHeader {
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..n_dip {
            hdr.push_str(&format!("\tS{s}"));
        }
        let lines = vec!["##fileformat=VCFv4.2".to_string(), hdr];
        VcfHeader::new_accept_all("src", &lines, &vec![true; n_dip])
    }

    fn rec(h: &VcfHeader, line: &str) -> Rc<dyn GTRec> {
        let p = VcfRecGTParser::new(h, line, &MarkerParser::new(true, true, true, true));
        Rc::new(BasicGTRec::from_parser(&p))
    }

    #[test]
    fn selects_common_low_missing_markers() {
        // 4 samples (8 haps). Two markers: m0 common (maf high), m1 rare (maf low).
        let h = header(4);
        let gt = BasicGT::new(vec![
            rec(&h, "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t0|1\t1|0\t1|0"), // 4 of 8 are alt
            rec(&h, "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|1"), // 1 of 8 alt
        ]);
        let mm = MarkerMap::create_min_dist(&PositionMap::new(1e-6), 0.0, gt.markers());
        // maf: m0 = 0.5 (>=0.1 -> usable), m1 = 0.125 (>=0.1 but only as computed)
        let maf = FloatArray::from_floats(&[0.5, 0.05]);
        let im = Ibs2Markers::new(&gt, &mm, &maf);
        assert_eq!(im.n_markers(), 2);
        // m0 usable (maf 0.5, no missing); m1 not usable (maf 0.05 < 0.1)
        assert_eq!(im.markers(0, 2), vec![0]);
    }

    #[test]
    fn missing_genotypes_exclude_marker() {
        let h = header(4);
        // marker with high maf but many missing genotypes -> excluded
        let gt = BasicGT::new(vec![rec(
            &h,
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t.|.\t.|.\t.|.",
        )]);
        let mm = MarkerMap::create_min_dist(&PositionMap::new(1e-6), 0.0, gt.markers());
        let maf = FloatArray::from_floats(&[0.5]);
        let im = Ibs2Markers::new(&gt, &mm, &maf);
        // 6 of 8 haps missing > maxMiss=ceil(0.1*8)=1 -> excluded
        assert_eq!(im.markers(0, 1), Vec::<i32>::new());
        let _ = im.step_starts().size();
    }
}
