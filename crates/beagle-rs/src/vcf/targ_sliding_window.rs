//! Port of `vcf/TargSlidingWindow.java` — a `SlidingWindow` over target-only VCF data.
//!
//! Java reads windows on a background thread into a 1-element blocking queue; the produced
//! window *sequence* is deterministic, so the Rust port reads them synchronously in
//! `next_window`.

use std::rc::Rc;

use crate::blbutil::{InputIt, SampleFileIt};
use crate::main_pkg::{Par, Pedigree};

use super::{
    genetic_map_from_file, marker_filter, sample_filter, to_lowmem_gt_rec, BasicGT, GTRec,
    GeneticMap, IntervalVcfIt, Marker, MarkerIndices, Samples, SlidingWindow, VcfIt, Window, GT,
};

type Rec = Box<dyn GTRec>;
type TargIt = Box<dyn SampleFileIt<Item = Rec>>;

/// Port of `vcf/TargSlidingWindow.java`.
pub struct TargSlidingWindow {
    targ_it: TargIt,
    ped: Pedigree,
    gen_map: Rc<dyn GeneticMap>,
    window_cm: f32,
    window_markers: i32,
    overlap_cm: f32,
    overlap_markers: i32,
    overlap: Vec<Rc<dyn GTRec>>,
    recs: Vec<Rc<dyn GTRec>>,
    next_rec: Option<Rc<dyn GTRec>>,
    window_index: i32,
    cum_targ_markers: i32,
    no_more_windows: bool,
}

fn to_rc(b: Box<dyn GTRec>) -> Rc<dyn GTRec> {
    Rc::from(b)
}

fn build_targ_it(par: &Par) -> TargIt {
    let n_buffered_blocks = par.nthreads() << 2;
    let it = InputIt::from_bgzip_file(par.gt(), n_buffered_blocks);
    let s_filter = sample_filter(par.excludesamples());
    let m_filter = marker_filter(par.excludemarkers());
    let vcf_it = VcfIt::create_filtered(it, s_filter, m_filter, to_lowmem_gt_rec);
    match par.chrom_int() {
        Some(ci) => Box::new(IntervalVcfIt::new(vcf_it, ci.clone())),
        None => Box::new(vcf_it),
    }
}

fn first_index_with_pos(targ_gt: &BasicGT, mut index: i32) -> i32 {
    let pos = targ_gt.marker(index).pos();
    while index > 0 && targ_gt.marker(index - 1).pos() == pos {
        index -= 1;
    }
    index
}

fn targ_overlap_start(
    gen_map: &dyn GeneticMap,
    overlap_cm: f32,
    overlap_markers: i32,
    targ_gt: &BasicGT,
    chrom_end: bool,
) -> i32 {
    if chrom_end {
        return targ_gt.n_markers();
    }
    let n_markers_m1 = targ_gt.n_markers() - 1;
    let marker: &Marker = targ_gt.marker(n_markers_m1);
    let end_gen_pos = gen_map.gen_pos_marker(marker);
    let start_gen_pos = end_gen_pos - overlap_cm as f64;
    let key = gen_map.base_pos(marker.chrom_index(), start_gen_pos);
    let mut low = (targ_gt.n_markers() - overlap_markers).max(0);
    let mut high = n_markers_m1;
    while low <= high {
        let mid = (low + high) >> 1;
        let mid_pos = targ_gt.marker(mid).pos();
        if mid_pos < key {
            low = mid + 1;
        } else if mid_pos > key {
            high = mid - 1;
        } else {
            return first_index_with_pos(targ_gt, mid);
        }
    }
    debug_assert!(high < low);
    first_index_with_pos(targ_gt, high.max(0))
}

impl TargSlidingWindow {
    /// `TargSlidingWindow.instance(Par par)`.
    pub fn instance(par: &Par) -> Self {
        assert!(par.ref_file().is_none(), "par.ref() is present");
        let mut targ_it = build_targ_it(par);
        let ped = Pedigree::new(targ_it.samples().clone(), par.ped());
        let gen_map: Rc<dyn GeneticMap> =
            Rc::from(genetic_map_from_file(par.map(), par.chrom_int()));
        let next_rec = targ_it.next().map(to_rc);
        assert!(next_rec.is_some(), "Error: no genotype data");
        let window_markers = par.window_markers();
        TargSlidingWindow {
            targ_it,
            ped,
            gen_map,
            window_cm: par.window(),
            window_markers,
            overlap_cm: par.overlap(),
            overlap_markers: window_markers >> 2,
            overlap: Vec::new(),
            recs: Vec::with_capacity(10_000),
            next_rec,
            window_index: 0,
            cum_targ_markers: 0,
            no_more_windows: false,
        }
    }

    fn next_end_cm(&self, next_marker: &Marker) -> f64 {
        let mut end_cm = self.gen_map.gen_pos_marker(next_marker);
        if self.overlap.is_empty() {
            end_cm += self.window_cm as f64;
        } else {
            end_cm += (self.window_cm - self.overlap_cm) as f64;
        }
        end_cm
    }

    fn read_window(&mut self, chrom_index: i32, end_pos: i32, window_index: i32) -> Window {
        let overlap_end = self.overlap.len() as i32;
        self.recs.clear();
        self.recs.append(&mut self.overlap);
        loop {
            let keep = match &self.next_rec {
                Some(rec) => {
                    let m = rec.marker();
                    m.chrom_index() == chrom_index
                        && m.pos() < end_pos
                        && (self.recs.len() as i32) < self.window_markers
                }
                None => false,
            };
            if !keep {
                break;
            }
            let rec = self.next_rec.take().expect("rec present");
            self.recs.push(rec);
            self.next_rec = self.targ_it.next().map(to_rc);
        }
        let targ_gt = BasicGT::new(self.recs.clone());
        let last_window = self.next_rec.is_none();
        let chrom_end = match &self.next_rec {
            None => true,
            Some(r) => r.marker().chrom_index() != chrom_index,
        };
        let overlap_start = targ_overlap_start(
            self.gen_map.as_ref(),
            self.overlap_cm,
            self.overlap_markers,
            &targ_gt,
            chrom_end,
        );
        let marker_indices =
            MarkerIndices::from_counts(overlap_end, overlap_start, targ_gt.n_markers());
        Window::new(
            self.gen_map.clone(),
            window_index,
            last_window,
            marker_indices,
            None,
            targ_gt,
        )
    }
}

impl SlidingWindow for TargSlidingWindow {
    fn targ_samples(&self) -> &Samples {
        self.targ_it.samples()
    }

    fn ped(&self) -> &Pedigree {
        &self.ped
    }

    fn gen_map(&self) -> &dyn GeneticMap {
        self.gen_map.as_ref()
    }

    fn cum_targ_markers(&self) -> i32 {
        self.cum_targ_markers
    }

    fn cum_markers(&self) -> i32 {
        self.cum_targ_markers
    }

    fn next_window(&mut self) -> Option<Window> {
        if self.no_more_windows {
            return None;
        }
        let (chrom_index, next_end_cm) = {
            let marker = self.next_rec.as_ref().expect("next_rec present").marker();
            (marker.chrom_index(), self.next_end_cm(marker))
        };
        let end_pos = self.gen_map.base_pos(chrom_index, next_end_cm);
        self.window_index += 1;
        let window = self.read_window(chrom_index, end_pos, self.window_index);
        let overlap_start = window.indices().overlap_start();
        self.overlap = self.recs[overlap_start as usize..].to_vec();
        if window.last_window() {
            self.no_more_windows = true;
        }
        let n_targ = window.indices().n_targ_markers();
        let targ_overlap_end = window.indices().targ_overlap_end();
        self.cum_targ_markers += n_targ - targ_overlap_end;
        Some(window)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_vcf(name: &str, body: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(name);
        std::fs::File::create(&path)
            .unwrap()
            .write_all(body.as_bytes())
            .unwrap();
        path
    }

    #[test]
    fn single_window_covers_close_markers() {
        let gt = write_vcf(
            "beagle_rs_tsw.vcf",
            "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\tS1\n\
chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0\n\
chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|0\t0|1\n\
chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t1|1\t0|0\n",
        );
        let par = Par::new(&[
            format!("gt={}", gt.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
            "seed=99999".to_string(),
        ]);
        let mut sw = TargSlidingWindow::instance(&par);
        assert_eq!(sw.targ_samples().size(), 2);
        assert_eq!(sw.ped().n_singles(), 2);

        let w = sw.next_window().expect("first window");
        assert_eq!(w.window_index(), 1);
        assert!(w.last_window());
        assert_eq!(w.targ_gt().n_markers(), 3);
        assert!(w.ref_gt().is_none());
        // all 3 target markers, no overlap with a (nonexistent) previous window
        assert_eq!(sw.cum_targ_markers(), 3);
        // no more windows
        assert!(sw.next_window().is_none());
    }
}
