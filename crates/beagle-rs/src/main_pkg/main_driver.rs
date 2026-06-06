//! Port of `main/Main.java` — the top-level Beagle driver: it parses parameters, opens the
//! sliding window over the input, and for each window phases (and optionally imputes) the
//! genotypes, writing the result through [`WindowWriter`] and statistics through [`RunStats`].
//!
//! Java spawns a producer thread for the sliding window and a `ForkJoinPool` for per-record
//! parallelism; both are sequential here (windows are produced synchronously, records are
//! compressed in index order), so `Locale.setDefault`, the `ForkJoinPool.parallelism` system
//! property, and `SlidingWindow.close()` have no Rust equivalent and are omitted. Output is
//! byte-identical because all parallel Java steps preserve encounter order.

use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::blbutil::{consts, Utilities};
use crate::imp::{ImpData, ImpLS};
use crate::jdk::Random;
use crate::phase::{FixedPhaseData, PhaseData, PhaseLS, Stage2Haps, SwapRate};
use crate::vcf::{RefTargSlidingWindow, SlidingWindow, TargSlidingWindow, Window, XRefGT, GT};

use super::{short_help, Par, RunStats, WindowWriter, COPYRIGHT, PROGRAM, VERSION};

use std::rc::Rc;

/// Port of `main/Main.java`.
struct Main {
    par: Par,
    sliding_wind: Box<dyn SlidingWindow>,
    run_stats: RunStats,
    window_writer: WindowWriter,
    /// generates a distinct seed for each window
    rand: Random,
}

/// `Main.main(String[] args)` — the Beagle program entry point.
pub fn main(args: &[String]) {
    // Java: Locale.setDefault(Locale.US) — the Rust port formats with US conventions directly.
    if args.is_empty() {
        println!("{PROGRAM} {VERSION}");
        println!("{COPYRIGHT}");
        println!("{}", Par::usage());
        std::process::exit(0);
    }
    let par = parameters(args);
    // Java sets ForkJoinPool.common.parallelism here; the Rust port is sequential.
    let mut run_stats = RunStats::new(&par);
    run_stats.print_start_info();

    let sliding_wind = sliding_window(&par);
    let window_writer = WindowWriter::new(&par, sliding_wind.targ_samples().clone());
    let mut main = Main::new(par, sliding_wind, run_stats, window_writer);
    main.phase_and_impute();
    let cum_targ = main.sliding_wind.cum_targ_markers();
    let cum_markers = main.sliding_wind.cum_markers();
    main.run_stats
        .print_summary_and_close(cum_targ, cum_markers);
    // Java closes winOut then slidingWind via try-with-resources (reverse declaration order).
    main.window_writer.close();
}

fn sliding_window(par: &Par) -> Box<dyn SlidingWindow> {
    if par.ref_file().is_none() {
        Box::new(TargSlidingWindow::instance(par))
    } else {
        Box::new(RefTargSlidingWindow::instance(par))
    }
}

/// `parameters(String[] args)` — builds `Par` and validates the output prefix and the
/// window/overlap relationship.
fn parameters(args: &[String]) -> Par {
    // warnings are printed in RunStats::print_start_info
    let par = Par::new(args);
    check_output_prefix(&par);
    if (par.window() as f64) < 1.1 * (par.overlap() as f64) {
        let s = format!(
            "{}{nl}{nl}ERROR: The \"window\" parameter must be at least 1.1 times the \"overlap\" parameter{nl}Exiting program.",
            short_help(),
            nl = consts::NL
        );
        Utilities::exit(&s);
    }
    par
}

fn check_output_prefix(par: &Par) {
    let out_prefix = Path::new(par.out());
    if out_prefix.is_dir() {
        let s = format!(
            "ERROR: \"out\" parameter cannot be a directory: \"{}\"",
            par.out()
        );
        Utilities::exit(&format!("{}{s}", Par::usage()));
    }

    let vcf_out = PathBuf::from(format!("{}.vcf.gz", par.out()));
    if let Some(r) = par.ref_file() {
        if vcf_out == r {
            let s = format!("ERROR: VCF output file equals input file: {}", r.display());
            Utilities::exit(&format!("{}{s}", Par::usage()));
        }
    }
    if vcf_out == par.gt() {
        let s = format!(
            "ERROR: VCF output file equals input file: {}",
            par.gt().display()
        );
        Utilities::exit(&format!("{}{s}", Par::usage()));
    }
}

impl Main {
    fn new(
        par: Par,
        sliding_wind: Box<dyn SlidingWindow>,
        run_stats: RunStats,
        window_writer: WindowWriter,
    ) -> Self {
        let rand = Random::new(par.seed());
        Main {
            par,
            sliding_wind,
            run_stats,
            window_writer,
            rand,
        }
    }

    /// `phaseAndImpute()`.
    fn phase_and_impute(&mut self) {
        let mut opt_window = self.sliding_wind.next_window();
        self.print_sample_summary(opt_window.as_ref());
        let mut overlap: Option<Rc<dyn GT>> = None;
        while let Some(window) = opt_window {
            let fpd = Rc::new(FixedPhaseData::new(
                &self.par,
                self.sliding_wind.ped(),
                &window,
                overlap.clone(),
            ));
            self.run_stats.print_window_update(&window, &fpd);
            let seed = self.rand.next_long();
            let pd = PhaseData::new(fpd.clone(), seed);
            overlap = if fpd.targ_gt().is_phased() {
                let phased_targ = Rc::new(XRefGT::from_phased_gt(
                    fpd.targ_gt().as_ref(),
                    self.par.nthreads(),
                ));
                self.print_window_xref(&window, phased_targ)
            } else {
                self.phase_stage1_variants(&pd);
                if fpd.stage1_targ_gt().n_markers() == fpd.targ_gt().n_markers() {
                    let phased_targ = Rc::new(pd.est_phase().phased_haps());
                    self.print_window_xref(&window, phased_targ)
                } else {
                    let stage2_haps = self.phase_stage2_variants(&pd);
                    self.print_window_stage2(&window, &stage2_haps)
                }
            };
            opt_window = self.sliding_wind.next_window();
        }
        // Java: slidingWind.close() — no Rust equivalent (synchronous producer).
    }

    fn print_sample_summary(&mut self, opt_window: Option<&Window>) {
        if let Some(window) = opt_window {
            let ped = self.sliding_wind.ped();
            self.run_stats.print_sample_summary(ped, window);
        }
    }

    /// `phaseStage1Variants(PhaseData pd)`.
    fn phase_stage1_variants(&mut self, pd: &PhaseData) {
        let n_its = self.par.burnin() + self.par.iterations();
        let max_burnin_swap_rate = 0.01_f64;
        while pd.it() < n_its {
            let t0 = Instant::now();
            PhaseLS::run_stage1(pd);
            self.run_stats
                .print_stage1_info(pd, t0.elapsed().as_nanos() as i64);
            pd.increment_it();
            let swap_rate = SwapRate::get_and_reset_swap_rate();
            if pd.it() < self.par.burnin() && swap_rate <= max_burnin_swap_rate {
                pd.advance_to_first_phasing_it();
            }
        }
    }

    /// `phaseStage2Variants(PhaseData pd)`.
    fn phase_stage2_variants(&mut self, pd: &PhaseData) -> Stage2Haps {
        let t0 = Instant::now();
        let stage2_haps = PhaseLS::run_stage2(pd);
        self.run_stats
            .print_stage2_info(t0.elapsed().as_nanos() as i64);
        stage2_haps
    }

    /// `printWindow(Window window, XRefGT phasedTarg)`.
    fn print_window_xref(
        &mut self,
        window: &Window,
        phased_targ: Rc<XRefGT>,
    ) -> Option<Rc<dyn GT>> {
        let impute = window.indices().n_markers() != window.indices().n_targ_markers();
        if !impute {
            let m_start = window.indices().prev_targ_splice();
            let m_end = window.indices().next_targ_splice();
            self.window_writer
                .print_phased_gt(phased_targ.as_ref(), m_start, m_end);
            self.phased_overlap(window, &phased_targ)
        } else {
            let t0 = Instant::now();
            let imp_data = ImpData::new(
                &self.par,
                window,
                phased_targ.clone() as Rc<dyn GT>,
                window.gen_map(),
            );
            let state_probs = ImpLS::state_probs(&imp_data);
            let m_start = window.indices().prev_splice();
            let m_end = window.indices().next_splice();
            self.window_writer
                .print_imputed(&imp_data, m_start, m_end, &state_probs);
            self.run_stats
                .imputation_nanos(t0.elapsed().as_nanos() as i64);
            self.run_stats.print_imputation_update();
            self.phased_overlap(window, &phased_targ)
        }
    }

    /// `printWindow(Window window, Stage2Haps stage2Haps)`.
    fn print_window_stage2(
        &mut self,
        window: &Window,
        stage2_haps: &Stage2Haps,
    ) -> Option<Rc<dyn GT>> {
        let impute = window.indices().n_markers() != window.indices().n_targ_markers();
        if !impute {
            let m_start = window.indices().prev_targ_splice();
            let m_end = window.indices().next_targ_splice();
            self.window_writer
                .print_phased_stage2(stage2_haps, m_start, m_end);
            let basic = stage2_haps.to_basic_gt(window.indices().targ_overlap_start(), m_end);
            Some(Rc::new(basic) as Rc<dyn GT>)
        } else {
            let t0 = Instant::now();
            let phased_hap_major = stage2_haps.to_basic_gt(0, window.targ_gt().n_markers());
            let phased_targ = Rc::new(XRefGT::from_phased_gt(
                &phased_hap_major,
                self.par.nthreads(),
            ));
            let imp_data = ImpData::new(
                &self.par,
                window,
                phased_targ.clone() as Rc<dyn GT>,
                window.gen_map(),
            );
            let state_probs = ImpLS::state_probs(&imp_data);
            let m_start = window.indices().prev_splice();
            let m_end = window.indices().next_splice();
            self.window_writer
                .print_imputed(&imp_data, m_start, m_end, &state_probs);
            self.run_stats
                .imputation_nanos(t0.elapsed().as_nanos() as i64);
            self.run_stats.print_imputation_update();
            self.phased_overlap(window, &phased_targ)
        }
    }

    /// `phasedOverlap(Window window, XRefGT phasedTarg)`.
    fn phased_overlap(&self, window: &Window, phased_targ: &XRefGT) -> Option<Rc<dyn GT>> {
        debug_assert!(phased_targ.is_phased());
        let mi = window.indices();
        let next_overlap = mi.targ_overlap_start();
        let next_splice = mi.next_targ_splice();
        let n_markers = next_splice - next_overlap;
        if n_markers == 0 {
            None
        } else {
            Some(Rc::new(phased_targ.restrict_to_xref(next_overlap, next_splice)) as Rc<dyn GT>)
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::{Read, Write};

    use flate2::read::MultiGzDecoder;

    fn write_temp(name: &str, body: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(name);
        std::fs::File::create(&path)
            .unwrap()
            .write_all(body.as_bytes())
            .unwrap();
        path
    }

    /// End-to-end smoke test: phase a small target-only (no reference) VCF and confirm the
    /// driver writes a valid BGZIP `<out>.vcf.gz` whose decompressed body is a phased VCF with
    /// one record per input marker and the correct sample columns.
    #[test]
    fn phases_target_only_vcf_end_to_end() {
        let gt = write_temp(
            "beagle_rs_main_gt.vcf",
            "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\tS1\tS2\tS3\n\
chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t1/0\t0/0\t1/1\n\
chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1\n\
chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\t0/0\t0/1\t1/0\n\
chr1\t400\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0\t0/1\n\
chr1\t500\t.\tT\tA\t.\tPASS\t.\tGT\t0/0\t0/1\t1/0\t1/1\n",
        );
        let out = std::env::temp_dir().join("beagle_rs_main_out");
        let args = vec![
            format!("gt={}", gt.display()),
            format!("out={}", out.display()),
            "nthreads=1".to_string(),
            "seed=99999".to_string(),
        ];
        super::main(&args);

        let vcf_gz = std::env::temp_dir().join("beagle_rs_main_out.vcf.gz");
        let bytes = std::fs::read(&vcf_gz).expect("output .vcf.gz exists");
        let mut decoded = String::new();
        MultiGzDecoder::new(&bytes[..])
            .read_to_string(&mut decoded)
            .expect("valid BGZIP");

        let header = decoded
            .lines()
            .find(|l| l.starts_with("#CHROM"))
            .expect("has #CHROM header");
        assert!(header.ends_with("FORMAT\tS0\tS1\tS2\tS3"));
        let data_lines: Vec<&str> = decoded
            .lines()
            .filter(|l| l.starts_with("chr1\t"))
            .collect();
        assert_eq!(data_lines.len(), 5, "one record per input marker");
        for line in &data_lines {
            // every output genotype must be phased ('|')
            let gt_fields = &line.split('\t').collect::<Vec<_>>()[9..];
            for g in gt_fields {
                assert!(g.contains('|'), "genotype not phased: {g}");
            }
        }
    }
}
