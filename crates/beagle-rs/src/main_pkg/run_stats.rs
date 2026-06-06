//! Port of `main/RunStats.java` — stores and prints analysis statistics to the `.log` file and
//! to standard output.
//!
//! Wall-clock fields (start/end time, elapsed times) are the documented `.log` normalization
//! exception; everything else (sample/marker counts, iteration labels, estimated parameters) is
//! byte-exact, including Java's `%,Nd` grouped integers, `%-Ns` left-justified labels, and the
//! `%7.1e` scientific format for the estimated error.

use std::io::Write;
use std::time::Instant;

use crate::blbutil::{consts, FileUtil, Utilities};
use crate::phase::{FixedPhaseData, PhaseData};
use crate::vcf::{Window, GT};

use super::{short_help, Par, Pedigree, PROGRAM};

/// `String.format("%,Nd", n)` — grouped integer (comma thousands separators), right-justified
/// to width `w`.
fn grouped_right(n: i64, w: usize) -> String {
    let neg = n < 0;
    let digits = n.unsigned_abs().to_string();
    let bytes = digits.as_bytes();
    let len = bytes.len();
    let mut g = String::with_capacity(len + len / 3 + 1);
    if neg {
        g.push('-');
    }
    for (i, b) in bytes.iter().enumerate() {
        if i > 0 && (len - i).is_multiple_of(3) {
            g.push(',');
        }
        g.push(*b as char);
    }
    format!("{g:>w$}")
}

/// `String.format("%7.1e", x)` — scientific, one fractional digit, signed ≥2-digit exponent,
/// right-justified to width 7.
fn sci_1e_w7(x: f32) -> String {
    let s = format!("{:.1e}", x as f64); // e.g. "1.0e-4", "0.0e0"
    let (mant, exp) = s.split_once('e').unwrap_or((s.as_str(), "0"));
    let exp_i: i32 = exp.parse().unwrap_or(0);
    let sign = if exp_i < 0 { '-' } else { '+' };
    let formatted = format!("{mant}e{sign}{:02}", exp_i.abs());
    format!("{formatted:>7}")
}

/// Port of `main/RunStats.java`.
pub struct RunStats {
    par: Par,
    log: Box<dyn Write>,
    start: Instant,
    total_phase_nanos: i64,
    impute_nanos: i64,
    total_impute_nanos: i64,
}

impl RunStats {
    /// `new RunStats(Par par)`.
    pub fn new(par: &Par) -> Self {
        let log = FileUtil::print_writer(std::path::Path::new(&format!("{}.log", par.out())));
        RunStats {
            par: par.clone(),
            log,
            start: Instant::now(),
            total_phase_nanos: 0,
            impute_nanos: 0,
            total_impute_nanos: 0,
        }
    }

    /// `par()`.
    pub fn par(&self) -> &Par {
        &self.par
    }

    /// `printStartInfo()`.
    pub fn print_start_info(&mut self) {
        let mut arg_list = self.par.args();
        if self.par.no_n_threads() {
            arg_list.push(format!("nthreads={}", self.par.nthreads()));
        }
        Utilities::duo_print(&mut *self.log, &format!("{}{}", short_help(), consts::NL));
        Utilities::duo_println(
            &mut *self.log,
            &format!("Start time: {}", Utilities::time_stamp()),
        );
        Utilities::duo_print(&mut *self.log, &Utilities::command_line(PROGRAM, &arg_list));
        if self.par.ped().is_some() {
            let s = format!(
                "{}WARNING: This version will not model duos or trios in the pedigree file",
                consts::NL
            );
            Utilities::duo_println(&mut *self.log, &s);
        }
        if self.par.map().is_none() {
            let s = format!(
                "{}No genetic map is specified: using 1 cM = 1 Mb",
                consts::NL
            );
            Utilities::duo_println(&mut *self.log, &s);
        }
        let _ = self.log.flush();
    }

    /// `printSampleSummary(Pedigree ped, Window window)`.
    pub fn print_sample_summary(&mut self, ped: &Pedigree, window: &Window) {
        let n_ref_samples = window.ref_gt().map(|r| r.n_samples()).unwrap_or(0);
        Utilities::duo_print(&mut *self.log, consts::NL);
        Utilities::duo_print(
            &mut *self.log,
            &format!(
                "Reference samples: {}\n",
                grouped_right(n_ref_samples as i64, 20)
            ),
        );
        Utilities::duo_print(
            &mut *self.log,
            &format!(
                "Study     samples: {}\n",
                grouped_right(window.targ_gt().n_samples() as i64, 20)
            ),
        );
        if self.par.ped().is_some() {
            Utilities::duo_print(&mut *self.log, "  ");
            Utilities::duo_print(&mut *self.log, &ped.n_singles().to_string());
            Utilities::duo_println(&mut *self.log, " singles");
            Utilities::duo_print(&mut *self.log, "  ");
            Utilities::duo_print(&mut *self.log, &ped.n_duos().to_string());
            Utilities::duo_println(&mut *self.log, " duos");
            Utilities::duo_print(&mut *self.log, "  ");
            Utilities::duo_print(&mut *self.log, &ped.n_trios().to_string());
            Utilities::duo_println(&mut *self.log, " trios");
        }
        let _ = self.log.flush();
    }

    /// `printWindowUpdate(Window window, FixedPhaseData fpd)`.
    pub fn print_window_update(&mut self, window: &Window, fpd: &FixedPhaseData) {
        let targ_markers = fpd.targ_gt().markers().clone();
        let stage1_targ_markers = fpd.stage1_targ_gt().markers().clone();
        let markers = window
            .ref_gt()
            .map(|r| r.markers().clone())
            .unwrap_or_else(|| targ_markers.clone());
        let first = markers.marker(0);
        let last = markers.marker(markers.size() - 1);
        let mut sb = String::with_capacity(30);
        sb.push_str(consts::NL);
        sb.push_str("Window ");
        sb.push_str(&fpd.window().to_string());
        sb.push_str(" [");
        let chr = first.chrom();
        if chr != consts::MISSING_DATA_STRING {
            sb.push_str(&chr);
            sb.push(consts::COLON);
        }
        sb.push_str(&first.pos().to_string());
        sb.push(consts::HYPHEN);
        if chr != last.chrom() {
            sb.push_str(&last.chrom());
            sb.push(consts::COLON);
        }
        sb.push_str(&last.pos().to_string());
        sb.push(']');
        sb.push_str(consts::NL);
        if window.ref_gt().is_some() {
            sb.push_str(&format!(
                "Reference markers: {}\n",
                grouped_right(markers.size() as i64, 20)
            ));
        }
        sb.push_str(&format!(
            "Study     markers: {}\n",
            grouped_right(targ_markers.size() as i64, 20)
        ));
        if stage1_targ_markers.size() != targ_markers.size() {
            sb.push_str(&format!(
                "Stage 1   markers: {}\n",
                grouped_right(stage1_targ_markers.size() as i64, 20)
            ));
        }
        Utilities::duo_print(&mut *self.log, &sb);
        let _ = self.log.flush();
    }

    /// `printSummaryAndClose(int nTargetMarkers, int nMarkers)`.
    pub fn print_summary_and_close(&mut self, n_target_markers: i32, n_markers: i32) {
        let total_time = self.start.elapsed().as_nanos() as i64;
        Utilities::duo_print(&mut *self.log, consts::NL);
        Utilities::duo_println(
            &mut *self.log,
            &format!("Cumulative Statistics:{}", consts::NL),
        );
        if n_target_markers != n_markers {
            Utilities::duo_print(
                &mut *self.log,
                &format!(
                    "Reference markers: {}\n",
                    grouped_right(n_markers as i64, 20)
                ),
            );
        }
        Utilities::duo_print(
            &mut *self.log,
            &format!(
                "Study     markers: {}\n\n",
                grouped_right(n_target_markers as i64, 20)
            ),
        );
        if self.total_phase_nanos > 1000 {
            self.duo_print_nanos("Haplotype phasing time:        ", self.total_phase_nanos);
        }
        if self.total_impute_nanos > 0 {
            self.duo_print_nanos("Imputation time:               ", self.total_impute_nanos);
        }
        self.duo_print_nanos("Total time:                    ", total_time);
        Utilities::duo_println(
            &mut *self.log,
            &format!("{}End time: {}", consts::NL, Utilities::time_stamp()),
        );
        Utilities::duo_println(&mut *self.log, &format!("{PROGRAM} finished"));
        let _ = self.log.flush();
    }

    /// `phaseNanos(long nanos)`.
    pub fn phase_nanos(&mut self, nanos: i64) {
        self.total_phase_nanos += nanos;
    }

    /// `imputationNanos(long nanos)`.
    pub fn imputation_nanos(&mut self, nanos: i64) {
        self.impute_nanos = nanos;
        self.total_impute_nanos += nanos;
    }

    /// `printImputationUpdate()`.
    pub fn print_imputation_update(&mut self) {
        Utilities::duo_print(&mut *self.log, consts::NL);
        let nanos = self.impute_nanos;
        self.duo_print_nanos("Imputation time:               ", nanos);
        let _ = self.log.flush();
    }

    /// `println(String msg)`.
    pub fn println(&mut self, msg: &str) {
        Utilities::duo_println(&mut *self.log, msg);
        let _ = self.log.flush();
    }

    /// `printStage1Info(PhaseData pd, long elapsedNanos)`.
    pub fn print_stage1_info(&mut self, pd: &PhaseData, elapsed_nanos: i64) {
        if pd.it() == self.par.burnin() && self.par.em() {
            self.print_estimated_parameters(pd.ne(), pd.p_mismatch());
        }
        self.phase_nanos(elapsed_nanos);
        let mut it = pd.it();
        let msg = if it < self.par.burnin() {
            if it == 0 {
                self.println("");
            }
            format!("Burnin  iteration {}:", it + 1)
        } else {
            it -= self.par.burnin();
            if it == 0 {
                self.println("");
            }
            format!("Phasing iteration {}:", it + 1)
        };
        self.duo_print_nanos(&format!("{msg:<31}"), elapsed_nanos);
    }

    /// `printStage2Info(long elapsedNanos)`.
    pub fn print_stage2_info(&mut self, elapsed_nanos: i64) {
        self.phase_nanos(elapsed_nanos);
        let msg = "Low frequency phasing:";
        self.duo_print_nanos(&format!("{msg:<31}"), elapsed_nanos);
    }

    /// `printEstimatedParameters(long ne, float pMismatch)`.
    pub fn print_estimated_parameters(&mut self, ne: i64, p_mismatch: f32) {
        Utilities::duo_println(&mut *self.log, "");
        Utilities::duo_print(&mut *self.log, &format!("{:<31}", "Estimated ne:"));
        Utilities::duo_println(&mut *self.log, &ne.to_string());
        Utilities::duo_print(&mut *self.log, &format!("{:<31}", "Estimated err:"));
        Utilities::duo_println(&mut *self.log, &sci_1e_w7(p_mismatch));
    }

    /// `duoPrintNanos(String message, long nanos)`.
    pub fn duo_print_nanos(&mut self, message: &str, nanos: i64) {
        Utilities::duo_print(&mut *self.log, message);
        Utilities::duo_println(&mut *self.log, &Utilities::elapsed_nanos(nanos));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grouped_and_scientific_formats() {
        // Java String.format("%,20d", n)
        assert_eq!(grouped_right(0, 20), "                   0");
        assert_eq!(grouped_right(1234, 20), "               1,234");
        assert_eq!(grouped_right(1234567, 20), "           1,234,567");
        // width narrower than the grouped value leaves it unpadded
        assert_eq!(grouped_right(1234567, 0), "1,234,567");
        assert_eq!(grouped_right(-1234, 0), "-1,234");
    }

    #[test]
    fn scientific_format_matches_java() {
        assert_eq!(sci_1e_w7(0.0001), "1.0e-04");
        assert_eq!(sci_1e_w7(0.0), "0.0e+00");
    }
}
