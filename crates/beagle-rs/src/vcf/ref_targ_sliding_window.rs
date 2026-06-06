//! Port of `vcf/RefTargSlidingWindow.java` — a `SlidingWindow` that merges reference and
//! target VCF records into overlapping cM windows. Reference-only markers are retained
//! (for imputation) interleaved with the shared target markers.

use std::rc::Rc;

use crate::blbutil::{consts, InputIt, SampleFileIt, Utilities};
use crate::bref::Bref3It;
use crate::main_pkg::{Par, Pedigree};

use super::{
    genetic_map_from_file, marker_filter, sample_filter, to_lowmem_gt_rec, BasicGT, GTRec,
    GeneticMap, IntervalVcfIt, Marker, MarkerIndices, RefGT, RefGTRec, RefIt, Samples,
    SlidingWindow, VcfIt, Window, GT,
};

const DEFAULT_BUFFER_SIZE: i32 = 1 << 10;

type TargIt = Box<dyn SampleFileIt<Item = Box<dyn GTRec>>>;
type RefItBox = Box<dyn SampleFileIt<Item = Box<dyn RefGTRec>>>;

/// Port of `vcf/RefTargSlidingWindow.java`.
pub struct RefTargSlidingWindow {
    targ_it: TargIt,
    ref_it: RefItBox,
    ped: Pedigree,
    gen_map: Rc<dyn GeneticMap>,
    window_cm: f32,
    window_markers: i32,
    overlap_cm: f32,
    overlap_markers: i32,
    impute: bool,
    targ_overlap: Vec<Rc<dyn GTRec>>,
    ref_overlap: Vec<Rc<dyn RefGTRec>>,
    in_targ_overlap: Vec<bool>,
    targ_recs: Vec<Rc<dyn GTRec>>,
    ref_recs: Vec<Rc<dyn RefGTRec>>,
    in_targ: Vec<bool>,
    next_targ_rec: Option<Rc<dyn GTRec>>,
    next_ref_rec: Option<Rc<dyn RefGTRec>>,
    window_index: i32,
    cum_targ_markers: i32,
    cum_ref_markers: i32,
    no_more_windows: bool,
}

fn to_rc(b: Box<dyn GTRec>) -> Rc<dyn GTRec> {
    Rc::from(b)
}

fn to_rc_ref(b: Box<dyn RefGTRec>) -> Rc<dyn RefGTRec> {
    Rc::from(b)
}

fn build_targ_it(par: &Par) -> TargIt {
    let n_buffered_blocks = par.nthreads() << 2;
    let it = InputIt::from_bgzip_file(par.gt(), n_buffered_blocks);
    let vcf_it = VcfIt::create_filtered(
        it,
        sample_filter(par.excludesamples()),
        marker_filter(par.excludemarkers()),
        to_lowmem_gt_rec,
    );
    match par.chrom_int() {
        Some(ci) => Box::new(IntervalVcfIt::new(vcf_it, ci.clone())),
        None => Box::new(vcf_it),
    }
}

fn build_ref_it(par: &Par) -> RefItBox {
    let ref_path = par.ref_file().expect("par.ref() present");
    let filename = ref_path.to_string_lossy();
    if filename.ends_with(".bref") {
        Utilities::exit(&format!(
            "{nl}ERROR: bref format (.bref) is not supported{nl}       \
             Reference files should be in bref3 format (.bref3)",
            nl = consts::NL
        ));
    }
    if filename.ends_with(".bref3") {
        let it = Bref3It::with_filters(
            Some(ref_path),
            sample_filter(par.excludesamples()),
            marker_filter(par.excludemarkers()),
        );
        match par.chrom_int() {
            Some(ci) => Box::new(IntervalVcfIt::new(it, ci.clone())),
            None => Box::new(it),
        }
    } else {
        if !filename.ends_with(".vcf")
            && !filename.ends_with(".vcf.gz")
            && !filename.ends_with(".vcf.bgz")
        {
            eprintln!(
                "{nl}ERROR: unrecognized reference filename extension: {nl}       \
                 expected \".bref3\", \".bref4\", \".vcf\", \".vcf.gz\", or \".vcf.bgz\"{nl}",
                nl = consts::NL
            );
        }
        let n_buffered_blocks = par.nthreads() << 2;
        let it = InputIt::from_bgzip_file(ref_path, n_buffered_blocks);
        let ref_it = RefIt::create_filtered(
            it,
            sample_filter(par.excludesamples()),
            marker_filter(par.excludemarkers()),
            DEFAULT_BUFFER_SIZE,
        );
        match par.chrom_int() {
            Some(ci) => Box::new(IntervalVcfIt::new(ref_it, ci.clone())),
            None => Box::new(ref_it),
        }
    }
}

fn first_index_with_pos(ref_gt: &RefGT, mut index: i32) -> i32 {
    let pos = ref_gt.marker(index).pos();
    while index > 0 && ref_gt.marker(index - 1).pos() == pos {
        index -= 1;
    }
    index
}

impl RefTargSlidingWindow {
    /// `RefTargSlidingWindow.instance(Par par)`.
    pub fn instance(par: &Par) -> Self {
        let mut targ_it = build_targ_it(par);
        let mut ref_it = build_ref_it(par);
        let ped = Pedigree::new(targ_it.samples().clone(), par.ped());
        let gen_map: Rc<dyn GeneticMap> =
            Rc::from(genetic_map_from_file(par.map(), par.chrom_int()));
        let next_targ_rec = targ_it.next().map(to_rc);
        let next_ref_rec = ref_it.next().map(to_rc_ref);
        assert!(
            next_targ_rec.is_some() && next_ref_rec.is_some(),
            "no genotype data"
        );
        let window_markers = par.window_markers();
        RefTargSlidingWindow {
            targ_it,
            ref_it,
            ped,
            gen_map,
            window_cm: par.window(),
            window_markers,
            overlap_cm: par.overlap(),
            overlap_markers: window_markers >> 2,
            impute: par.impute(),
            targ_overlap: Vec::new(),
            ref_overlap: Vec::new(),
            in_targ_overlap: Vec::new(),
            targ_recs: Vec::with_capacity(10_000),
            ref_recs: Vec::with_capacity(10_000),
            in_targ: Vec::with_capacity(10_000),
            next_targ_rec,
            next_ref_rec,
            window_index: 0,
            cum_targ_markers: 0,
            cum_ref_markers: 0,
            no_more_windows: false,
        }
    }

    fn advance_ref_it_to_chrom(&mut self, chrom_index: i32) {
        loop {
            match &self.next_ref_rec {
                Some(rec) if rec.marker().chrom_index() != chrom_index => {
                    match self.ref_it.next() {
                        Some(r) => self.next_ref_rec = Some(to_rc_ref(r)),
                        None => break, // refIt exhausted; keep last record
                    }
                }
                _ => break,
            }
        }
    }

    fn next_end_cm(&self, next_ref_marker: &Marker) -> f64 {
        let mut end_cm = self.gen_map.gen_pos_marker(next_ref_marker);
        if self.ref_overlap.is_empty() {
            end_cm += self.window_cm as f64;
        } else {
            end_cm += (self.window_cm - self.overlap_cm) as f64;
        }
        end_cm
    }

    fn reset_lists(&mut self) {
        self.targ_recs.clear();
        self.ref_recs.clear();
        self.in_targ.clear();
        self.targ_recs.append(&mut self.targ_overlap);
        self.ref_recs.append(&mut self.ref_overlap);
        self.in_targ.append(&mut self.in_targ_overlap);
    }

    fn read_window(&mut self, chrom_index: i32, end_pos: i32, window_index: i32) -> Window {
        let ref_overlap_end = self.ref_overlap.len() as i32;
        self.reset_lists();
        // not a `while let`: the body mutates `self.next_targ_rec`/`next_ref_rec`, so the
        // match's borrow must end before each iteration's mutations.
        #[allow(clippy::while_let_loop)]
        loop {
            let (targ_marker, targ_pos) = match &self.next_targ_rec {
                Some(rec) => {
                    let m = rec.marker();
                    if !(m.chrom_index() == chrom_index
                        && m.pos() < end_pos
                        && (self.ref_recs.len() as i32) < self.window_markers)
                    {
                        break;
                    }
                    (m.clone(), m.pos())
                }
                None => break,
            };
            loop {
                let advance = match &self.next_ref_rec {
                    Some(r) => {
                        let rm = r.marker();
                        rm.chrom_index() == chrom_index
                            && (rm.pos() < targ_pos || (rm.pos() == targ_pos && targ_marker != *rm))
                    }
                    None => false,
                };
                if !advance {
                    break;
                }
                if self.impute {
                    self.ref_recs
                        .push(self.next_ref_rec.clone().expect("ref rec"));
                    self.in_targ.push(false);
                }
                self.next_ref_rec = self.ref_it.next().map(to_rc_ref);
            }
            let matches = matches!(&self.next_ref_rec, Some(r) if *r.marker() == targ_marker);
            if matches {
                self.targ_recs
                    .push(self.next_targ_rec.clone().expect("targ rec"));
                self.ref_recs
                    .push(self.next_ref_rec.clone().expect("ref rec"));
                self.in_targ.push(true);
                self.next_ref_rec = self.ref_it.next().map(to_rc_ref);
            }
            self.next_targ_rec = self.targ_it.next().map(to_rc);
        }
        if self.impute {
            loop {
                let keep = match &self.next_ref_rec {
                    Some(r) => {
                        let m = r.marker();
                        m.chrom_index() == chrom_index
                            && m.pos() < end_pos
                            && (self.ref_recs.len() as i32) < self.window_markers
                    }
                    None => false,
                };
                if !keep {
                    break;
                }
                self.ref_recs
                    .push(self.next_ref_rec.clone().expect("ref rec"));
                self.in_targ.push(false);
                self.next_ref_rec = self.ref_it.next().map(to_rc_ref);
            }
        }
        self.build_window(ref_overlap_end, window_index, chrom_index, end_pos)
    }

    fn build_window(
        &self,
        ref_overlap_end: i32,
        window_index: i32,
        chrom_index: i32,
        end_pos: i32,
    ) -> Window {
        assert!(
            !self.targ_recs.is_empty() && !self.ref_recs.is_empty(),
            "{}",
            self.empty_window_error_message(chrom_index, end_pos)
        );
        let ref_gt = RefGT::from_recs(self.ref_recs.clone());
        let targ_gt = BasicGT::new(self.targ_recs.clone());
        let last_window = self.next_targ_rec.is_none() || self.next_ref_rec.is_none();
        let marker_indices = self.marker_indices(&ref_gt, last_window, ref_overlap_end, end_pos);
        Window::new(
            self.gen_map.clone(),
            window_index,
            last_window,
            marker_indices,
            Some(ref_gt),
            targ_gt,
        )
    }

    fn empty_window_error_message(&self, chrom_index: i32, end_pos: i32) -> String {
        use crate::beagleutil::ChromIds;
        if self.ref_recs.is_empty() {
            format!(
                "The window ending at {}:{end_pos}{nl}contains no reference markers{nl}\
                 Do the reference and target VCF files contain the same{nl}\
                 chromosomes in the same order?{nl}",
                ChromIds::instance().id(chrom_index),
                nl = consts::NL
            )
        } else {
            format!(
                "The reference and target VCF files contain no markers in common in the window: \
                 {nl}{}:{}-{end_pos}{nl}Do both VCF files share any markers in this window?{nl}\
                 Do both VCF files contain the same chromosomes in the same order?{nl}",
                ChromIds::instance().id(chrom_index),
                self.ref_recs[0].marker().pos(),
                nl = consts::NL
            )
        }
    }

    fn marker_indices(
        &self,
        ref_gt: &RefGT,
        last_window: bool,
        ref_overlap_end: i32,
        end_pos: i32,
    ) -> MarkerIndices {
        let chrom_end = last_window
            || (self.ref_recs[0].marker().chrom_index()
                != self
                    .next_ref_rec
                    .as_ref()
                    .expect("ref rec")
                    .marker()
                    .chrom_index());
        let ref_overlap_start = self.overlap_start(ref_gt, chrom_end, end_pos);
        MarkerIndices::from_in_targ(&self.in_targ, ref_overlap_end, ref_overlap_start)
    }

    fn overlap_start(&self, ref_gt: &RefGT, chrom_end: bool, end_pos: i32) -> i32 {
        if chrom_end {
            return ref_gt.n_markers();
        }
        let n_markers_m1 = ref_gt.n_markers() - 1;
        let chrom_index = ref_gt.marker(n_markers_m1).chrom_index();
        let end_gen_pos = self.gen_map.gen_pos(chrom_index, end_pos - 1);
        let start_gen_pos = end_gen_pos - self.overlap_cm as f64;
        let key = self.gen_map.base_pos(chrom_index, start_gen_pos);
        let mut low = (ref_gt.n_markers() - self.overlap_markers).max(0);
        let mut high = n_markers_m1;
        while low <= high {
            let mid = (low + high) >> 1;
            let mid_pos = ref_gt.marker(mid).pos();
            if mid_pos < key {
                low = mid + 1;
            } else if mid_pos > key {
                high = mid - 1;
            } else {
                return first_index_with_pos(ref_gt, mid);
            }
        }
        debug_assert!(high < low);
        first_index_with_pos(ref_gt, high.max(0))
    }
}

impl SlidingWindow for RefTargSlidingWindow {
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
        self.cum_ref_markers
    }

    fn next_window(&mut self) -> Option<Window> {
        if self.no_more_windows {
            return None;
        }
        let chrom_index = self
            .next_targ_rec
            .as_ref()
            .expect("next_targ_rec present")
            .marker()
            .chrom_index();
        self.advance_ref_it_to_chrom(chrom_index);
        let end_cm = {
            let m = self
                .next_ref_rec
                .as_ref()
                .expect("next_ref_rec present")
                .marker();
            self.next_end_cm(m)
        };
        let end_pos = self.gen_map.base_pos(chrom_index, end_cm);
        self.window_index += 1;
        let window = self.read_window(chrom_index, end_pos, self.window_index);

        let targ_overlap_start = window.indices().targ_overlap_start() as usize;
        let ref_overlap_start = window.indices().overlap_start() as usize;
        self.targ_overlap = self.targ_recs[targ_overlap_start..].to_vec();
        self.ref_overlap = self.ref_recs[ref_overlap_start..].to_vec();
        self.in_targ_overlap = self.in_targ[ref_overlap_start..].to_vec();

        if window.last_window() {
            self.no_more_windows = true;
        }
        let indices = window.indices();
        self.cum_targ_markers += indices.n_targ_markers() - indices.targ_overlap_end();
        self.cum_ref_markers += indices.n_markers() - indices.overlap_end();
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
    fn merges_ref_and_target_markers() {
        // reference has markers 100,150,200; target shares 100,200 (150 is ref-only).
        let ref_vcf = write_vcf(
            "beagle_rs_rtsw_ref.vcf",
            "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tR0\tR1\n\
chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0\n\
chr1\t150\t.\tG\tT\t.\tPASS\t.\tGT\t0|0\t1|1\n\
chr1\t200\t.\tA\tG\t.\tPASS\t.\tGT\t1|1\t0|0\n",
        );
        let gt_vcf = write_vcf(
            "beagle_rs_rtsw_gt.vcf",
            "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\n\
chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\n\
chr1\t200\t.\tA\tG\t.\tPASS\t.\tGT\t1|0\n",
        );
        let par = Par::new(&[
            format!("gt={}", gt_vcf.display()),
            format!("ref={}", ref_vcf.display()),
            "out=o".to_string(),
            "nthreads=1".to_string(),
            "seed=99999".to_string(),
        ]);
        let mut sw = RefTargSlidingWindow::instance(&par);
        assert_eq!(sw.targ_samples().size(), 1);

        let w = sw.next_window().expect("window");
        assert!(w.last_window());
        assert!(w.ref_gt().is_some());
        assert_eq!(w.ref_gt().unwrap().n_markers(), 3); // 100,150,200
        assert_eq!(w.targ_gt().n_markers(), 2); // 100,200
                                                // marker->target map: 100->0, 150->-1 (ref-only), 200->1
        assert_eq!(w.indices().marker_to_targ_marker(), vec![0, -1, 1]);
        assert_eq!(sw.cum_markers(), 3); // ref markers
        assert_eq!(sw.cum_targ_markers(), 2);
        assert!(sw.next_window().is_none());
    }
}
