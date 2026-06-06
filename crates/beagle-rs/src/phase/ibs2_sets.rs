//! Port of `phase/Ibs2Sets.java` — partitions markers into steps and stores, per step, the
//! sets of target samples whose genotypes are consistent with IBS2 (identical unordered
//! genotypes, not all homozygous).

use crate::ints::IntList;
use crate::vcf::GT;

use super::{Ibs2Markers, SampleSeg};

const MAX_MISS_STEP_FREQ: f32 = 0.1;

/// Port of `phase/Ibs2Sets.java`.
pub struct Ibs2Sets {
    n_targ_samples: i32,
    n_markers_m1: i32,
    window_starts: Vec<i32>,
    ibs2_sets: Vec<Vec<Vec<i32>>>, // [window][targ_sample][ibs2_samples]
}

struct SampClust {
    samples: Vec<i32>,
    is_homozygous: bool,
}

fn get_gt(m: i32, s: i32, targ_gt: &dyn GT) -> i32 {
    let hap1 = s << 1;
    let a1 = targ_gt.allele(m, hap1);
    let a2 = targ_gt.allele(m, hap1 | 0b1);
    if a1 < 0 || a2 < 0 {
        -1
    } else if a1 <= a2 {
        ((a2 * (a2 + 1)) >> 1) + a1
    } else {
        ((a1 * (a1 + 1)) >> 1) + a2
    }
}

fn is_hom(prev_is_hom: bool, n_alleles: i32) -> Vec<bool> {
    let n_gt = ((n_alleles * (n_alleles + 1)) >> 1) as usize;
    let mut is_hom = vec![false; n_gt];
    if prev_is_hom {
        for a in 0..n_alleles {
            is_hom[(((a * (a + 1)) >> 1) + a) as usize] = true;
        }
    }
    is_hom
}

fn init_cluster(targ_gt: &dyn GT, step_markers: &[i32]) -> SampClust {
    let n_targ_samples = targ_gt.n_samples();
    let mut miss_cnt = vec![0i32; n_targ_samples as usize];
    for &m in step_markers {
        for s in 0..n_targ_samples {
            let hap1 = s << 1;
            if targ_gt.allele(m, hap1) == -1 || targ_gt.allele(m, hap1 | 0b1) == -1 {
                miss_cnt[s as usize] += 1;
            }
        }
    }
    let max_miss = ((MAX_MISS_STEP_FREQ * step_markers.len() as f32) as f64).floor() as i32;
    let init: Vec<i32> = (0..n_targ_samples)
        .filter(|&s| miss_cnt[s as usize] <= max_miss)
        .collect();
    SampClust {
        samples: init,
        is_homozygous: true,
    }
}

fn partition(targ_gt: &dyn GT, parent: &SampClust, m: i32) -> Vec<SampClust> {
    let n_alleles = targ_gt.marker(m).n_alleles();
    let n_gt = ((n_alleles * (n_alleles + 1)) >> 1) as usize;
    let mut gt_to_list: Vec<Option<IntList>> = (0..n_gt).map(|_| None).collect();
    let is_hom = is_hom(parent.is_homozygous, n_alleles);
    let mut missing = IntList::with_capacity(32);
    for &s in &parent.samples {
        let gt_index = get_gt(m, s, targ_gt);
        if gt_index < 0 {
            missing.add(s);
            for l in gt_to_list.iter_mut().flatten() {
                l.add(s);
            }
        } else {
            let gi = gt_index as usize;
            if gt_to_list[gi].is_none() {
                let mut l = IntList::new();
                for j in 0..missing.size() {
                    l.add(missing.get(j));
                }
                gt_to_list[gi] = Some(l);
            }
            gt_to_list[gi].as_mut().expect("gt list").add(s);
        }
    }
    let mut out = Vec::new();
    for (i, slot) in gt_to_list.into_iter().enumerate() {
        if let Some(l) = slot {
            if l.size() > 1 {
                out.push(SampClust {
                    samples: l.to_array(),
                    is_homozygous: is_hom[i],
                });
            }
        }
    }
    out
}

fn results(ibd2_lists: &[SampClust], n_targ_samples: i32) -> Vec<Vec<i32>> {
    let mut results: Vec<Vec<i32>> = vec![Vec::new(); n_targ_samples as usize];
    for cluster in ibd2_lists {
        if !cluster.is_homozygous {
            let ia = &cluster.samples;
            for &s in ia {
                if results[s as usize].is_empty() {
                    results[s as usize] = ia.clone();
                } else {
                    // sample can be in >1 list due to missing genotypes
                    let mut merged: Vec<i32> = results[s as usize]
                        .iter()
                        .chain(ia.iter())
                        .copied()
                        .collect();
                    merged.sort_unstable();
                    merged.dedup();
                    results[s as usize] = merged;
                }
            }
        }
    }
    results
}

fn ibs2_sets_for_window(
    targ_gt: &dyn GT,
    ibs2_markers: &Ibs2Markers,
    w_starts: &[i32],
    w: usize,
) -> Vec<Vec<i32>> {
    let start = w_starts[w];
    let end = if w + 1 < w_starts.len() {
        w_starts[w + 1]
    } else {
        targ_gt.n_markers()
    };
    let step_markers = ibs2_markers.markers(start, end);
    let mut partition_v = vec![init_cluster(targ_gt, &step_markers)];
    for &m in &step_markers {
        let mut next = Vec::new();
        for sc in &partition_v {
            next.extend(partition(targ_gt, sc, m));
        }
        partition_v = next;
    }
    results(&partition_v, targ_gt.n_samples())
}

impl Ibs2Sets {
    /// `new Ibs2Sets(GT targGT, Ibs2Markers ibs2Markers)`.
    pub fn new(targ_gt: &dyn GT, ibs2_markers: &Ibs2Markers) -> Self {
        assert!(
            targ_gt.n_markers() == ibs2_markers.n_markers(),
            "{}",
            ibs2_markers.n_markers()
        );
        let window_starts = ibs2_markers.step_starts().to_array();
        let ibs2_sets: Vec<Vec<Vec<i32>>> = (0..window_starts.len())
            .map(|w| ibs2_sets_for_window(targ_gt, ibs2_markers, &window_starts, w))
            .collect();
        Ibs2Sets {
            n_targ_samples: targ_gt.n_samples(),
            n_markers_m1: targ_gt.n_markers() - 1,
            window_starts,
            ibs2_sets,
        }
    }

    /// `nTargSamples()`.
    pub fn n_targ_samples(&self) -> i32 {
        self.n_targ_samples
    }

    /// `segList(int sample)` — the IBS2 segments for the sample.
    pub fn seg_list(&self, sample: i32) -> Vec<SampleSeg> {
        let mut list = Vec::new();
        for (w, sets) in self.ibs2_sets.iter().enumerate() {
            let ia = &sets[sample as usize];
            if !ia.is_empty() {
                let start = self.window_starts[w];
                let incl_end = if w + 1 < self.window_starts.len() {
                    self.window_starts[w + 1] - 1
                } else {
                    self.n_markers_m1
                };
                for &s2 in ia {
                    if s2 != sample {
                        list.push(SampleSeg::new(s2, start, incl_end));
                    }
                }
            }
        }
        list
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blbutil::FloatArray;
    use crate::vcf::{
        BasicGT, BasicGTRec, GTRec, MarkerMap, MarkerParser, PositionMap, VcfHeader,
        VcfRecGTParser, HEADER_PREFIX,
    };
    use std::rc::Rc;

    #[test]
    fn detects_ibs2_between_identical_het_samples() {
        // 3 samples, 60 het markers spaced 50 kb apart. Samples 0 & 1 are het (0|1) at
        // every marker (IBS2 with each other); sample 2 is homozygous (0|0).
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
        let im = Ibs2Markers::new(&gt, &map, &maf);
        let sets = Ibs2Sets::new(&gt, &im);

        assert_eq!(sets.n_targ_samples(), 3);
        // sample 0 is IBS2 with sample 1
        let seg0 = sets.seg_list(0);
        assert!(seg0.iter().any(|ss| ss.sample() == 1));
        // sample 2 (homozygous) is in no IBS2 set
        assert!(sets.seg_list(2).is_empty());
    }
}
