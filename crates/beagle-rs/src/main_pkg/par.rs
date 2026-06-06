//! Port of `main/Par.java` — the immutable parameter set parsed from command-line args.

use std::path::{Path, PathBuf};

use crate::beagleutil::ChromInterval;
use crate::blbutil::{consts, validate};

// default phasing parameters
const D_BURNIN: i32 = 3;
const D_ITERATIONS: i32 = 12;
const D_INITIAL_LR: f32 = 100_000.0;
const D_PHASE_STATES: i32 = 280;
const D_STEP_SCALE: f32 = 3.0;
const D_RARE: f32 = 0.002;

// default imputation parameters
const D_IMPUTE: bool = true;
const D_IMP_STATES: i32 = 1600;
const D_IMP_SEGMENT: f32 = 6.0;
const D_IMP_STEP: f32 = 0.1;
const D_IMP_NSTEPS: i32 = 7;
const D_CLUSTER: f32 = 0.005;
const D_AP: bool = false;
const D_GP: bool = false;

// default general parameters
const D_EM: bool = true;
const D_NE: i32 = 100_000;
const D_WINDOW: f32 = 40.0;
const D_WINDOW_MARKERS: i32 = 4_000_000;
const D_OVERLAP: f32 = 2.0;
const D_SEED: i32 = -99999;
const D_NTHREADS: i32 = i32::MAX;
const D_BUFFER: f32 = 1.0;

/// Port of `main/Par.java`.
#[derive(Clone)]
pub struct Par {
    args: Vec<String>,
    no_n_threads: bool,

    gt: PathBuf,
    ref_file: Option<PathBuf>,
    out: String,
    map: Option<PathBuf>,
    chrom_int: Option<ChromInterval>,
    excludesamples: Option<PathBuf>,
    excludemarkers: Option<PathBuf>,

    burnin: i32,
    iterations: i32,
    initial_lr: f32,
    phase_states: i32,
    step_scale: f32,
    rare: f32,

    impute: bool,
    imp_states: i32,
    imp_segment: f32,
    imp_step: f32,
    imp_nsteps: i32,
    cluster: f32,
    ap: bool,
    gp: bool,

    em: bool,
    ne: f32,
    err: f32,
    window: f32,
    window_markers: i32,
    overlap: f32,
    seed: i64,
    nthreads: i32,
    buffer: f32,

    truth: Option<PathBuf>,
}

fn parse_chrom_int(s: Option<&str>) -> Option<ChromInterval> {
    let ci = s.and_then(ChromInterval::parse);
    if let Some(str) = s {
        if !str.is_empty() && ci.is_none() {
            panic!("Invalid chrom parameter: {str}");
        }
    }
    ci
}

impl Par {
    /// `new Par(String[] args)`.
    pub fn new(args: &[String]) -> Self {
        let imax = i32::MAX;
        let lmin = i64::MIN;
        let lmax = i64::MAX;
        // Java `Float.MIN_VALUE` is the smallest positive *subnormal* (2^-149), not Rust's
        // `f32::MIN_POSITIVE` (smallest normal). Reproduce it exactly via the bit pattern.
        let fmin = f32::from_bits(1);
        let fmax = f32::MAX;

        let mut m = validate::args_to_map(args, '=');

        // data input/output parameters
        let gt =
            validate::get_file(validate::string_arg("gt", &mut m, true, None, None).as_deref())
                .expect("required gt file");
        let ref_file =
            validate::get_file(validate::string_arg("ref", &mut m, false, None, None).as_deref());
        let out = validate::string_arg("out", &mut m, true, None, None).expect("required out");
        // ped is read (and discarded) to consume the arg; Par.ped() always returns null
        let _ped =
            validate::get_file(validate::string_arg("ped", &mut m, false, None, None).as_deref());
        let map =
            validate::get_file(validate::string_arg("map", &mut m, false, None, None).as_deref());
        let chrom_int =
            parse_chrom_int(validate::string_arg("chrom", &mut m, false, None, None).as_deref());
        let excludesamples = validate::get_file(
            validate::string_arg("excludesamples", &mut m, false, None, None).as_deref(),
        );
        let excludemarkers = validate::get_file(
            validate::string_arg("excludemarkers", &mut m, false, None, None).as_deref(),
        );

        // phasing parameters
        let burnin = validate::int_arg("burnin", &mut m, false, D_BURNIN, 1, imax);
        let iterations = validate::int_arg("iterations", &mut m, false, D_ITERATIONS, 1, imax);
        let initial_lr = validate::float_arg("initial-lr", &mut m, false, D_INITIAL_LR, 1.0, fmax);
        let phase_states =
            validate::int_arg("phase-states", &mut m, false, D_PHASE_STATES, 1, imax);
        let step_scale = validate::float_arg("step-scale", &mut m, false, D_STEP_SCALE, fmin, fmax);
        let rare = validate::float_arg("rare", &mut m, false, D_RARE, fmin, 0.5);

        // imputation parameters
        let impute = validate::boolean_arg("impute", &mut m, false, D_IMPUTE);
        let imp_states = validate::int_arg("imp-states", &mut m, false, D_IMP_STATES, 1, imax);
        let imp_segment =
            validate::float_arg("imp-segment", &mut m, false, D_IMP_SEGMENT, fmin, fmax);
        let imp_step = validate::float_arg("imp-step", &mut m, false, D_IMP_STEP, fmin, fmax);
        let imp_nsteps = validate::int_arg("imp-nsteps", &mut m, false, D_IMP_NSTEPS, 1, imax);
        let cluster = validate::float_arg("cluster", &mut m, false, D_CLUSTER, 0.0, fmax);
        let ap = validate::boolean_arg("ap", &mut m, false, D_AP);
        let gp = validate::boolean_arg("gp", &mut m, false, D_GP);

        // general parameters
        let em = validate::boolean_arg("em", &mut m, false, D_EM);
        let ne = validate::float_arg("ne", &mut m, false, D_NE as f32, fmin, fmax);
        let err = validate::float_arg("err", &mut m, false, -fmin, -fmin, fmax);
        let window = validate::float_arg("window", &mut m, false, D_WINDOW, fmin, fmax);
        let window_markers = validate::int_arg(
            "window-markers",
            &mut m,
            false,
            D_WINDOW_MARKERS,
            100_000,
            imax,
        );
        let overlap = validate::float_arg("overlap", &mut m, false, D_OVERLAP, fmin, fmax);
        let buffer = validate::float_arg("buffer", &mut m, false, D_BUFFER, fmin, fmax);
        let seed = validate::long_arg("seed", &mut m, false, D_SEED as i64, lmin, lmax);
        let raw_n_threads = validate::int_arg("nthreads", &mut m, false, D_NTHREADS, 1, imax);
        let no_n_threads = raw_n_threads == D_NTHREADS;
        let nthreads = if no_n_threads {
            available_processors()
        } else {
            raw_n_threads
        };

        // undocumented parameters
        let truth =
            validate::get_file(validate::string_arg("truth", &mut m, false, None, None).as_deref());

        validate::confirm_empty_map(&m);

        Par {
            args: args.to_vec(),
            no_n_threads,
            gt,
            ref_file,
            out,
            map,
            chrom_int,
            excludesamples,
            excludemarkers,
            burnin,
            iterations,
            initial_lr,
            phase_states,
            step_scale,
            rare,
            impute,
            imp_states,
            imp_segment,
            imp_step,
            imp_nsteps,
            cluster,
            ap,
            gp,
            em,
            ne,
            err,
            window,
            window_markers,
            overlap,
            seed,
            nthreads,
            buffer,
            truth,
        }
    }

    /// `args()`.
    pub fn args(&self) -> Vec<String> {
        self.args.clone()
    }

    /// `noNThreads()`.
    pub fn no_n_threads(&self) -> bool {
        self.no_n_threads
    }

    /// `usage()` — a description of the Beagle command-line arguments.
    pub fn usage() -> String {
        let nl = consts::NL;
        let cmd = super::COMMAND;
        let window = java_float_str(D_WINDOW);
        let overlap = java_float_str(D_OVERLAP);
        let cluster = java_float_str(D_CLUSTER);
        // Built as a line array joined by `nl` (a `\`-continuation in a single literal would eat
        // the leading two-space indent of each argument line).
        let lines: [String; 32] = [
            format!("Usage: {cmd} [arguments]"),
            String::new(),
            "data parameters ...".to_string(),
            "  gt=<VCF file with GT FORMAT field>                 (required)".to_string(),
            "  ref=<bref3 or VCF file with phased genotypes>      (optional)".to_string(),
            "  out=<output file prefix>                           (required)".to_string(),
            "  map=<PLINK map file with cM units>                 (optional)".to_string(),
            "  chrom=<[chrom] or [chrom]:[start]-[end]>           (optional)".to_string(),
            "  excludesamples=<file with 1 sample ID per line>    (optional)".to_string(),
            "  excludemarkers=<file with 1 marker ID per line>    (optional)".to_string(),
            String::new(),
            "phasing parameters ...".to_string(),
            format!("  burnin=<max burnin iterations>                     (default={D_BURNIN})"),
            format!(
                "  iterations=<phasing iterations>                    (default={D_ITERATIONS})"
            ),
            format!(
                "  phase-states=<model states for phasing>            (default={D_PHASE_STATES})"
            ),
            String::new(),
            "imputation parameters ...".to_string(),
            format!("  impute=<impute ungenotyped markers (true/false)>   (default={D_IMPUTE})"),
            format!(
                "  imp-states=<model states for imputation>           (default={D_IMP_STATES})"
            ),
            format!("  cluster=<max cM in a marker cluster>               (default={cluster})"),
            format!("  ap=<print posterior allele probabilities>          (default={D_AP})"),
            format!("  gp=<print posterior genotype probabilities>        (default={D_GP})"),
            String::new(),
            "general parameters ...".to_string(),
            format!("  ne=<effective population size>                     (default={D_NE})"),
            "  err=<allele mismatch probability>                  (default: data dependent)"
                .to_string(),
            format!("  em=<estimate ne and err parameters (true/false)>   (default={D_EM})"),
            format!("  window=<window length in cM>                       (default={window})"),
            format!(
                "  window-markers=<maximum markers per window>        (default={D_WINDOW_MARKERS})"
            ),
            format!("  overlap=<window overlap in cM>                     (default={overlap})"),
            format!("  seed=<random seed>                                 (default={D_SEED})"),
            "  nthreads=<number of threads>                       (default: machine dependent)"
                .to_string(),
        ];
        // Java terminates the last line and adds a trailing blank line (`... + nl + nl`).
        format!("{}{nl}{nl}", lines.join(nl))
    }

    // data parameters
    /// `gt()`.
    pub fn gt(&self) -> &Path {
        &self.gt
    }
    /// `ref()`.
    pub fn ref_file(&self) -> Option<&Path> {
        self.ref_file.as_deref()
    }
    /// `out()`.
    pub fn out(&self) -> &str {
        &self.out
    }
    /// `ped()` — always `None` (Java `Par.ped()` returns `null`).
    pub fn ped(&self) -> Option<&Path> {
        None
    }
    /// `map()`.
    pub fn map(&self) -> Option<&Path> {
        self.map.as_deref()
    }
    /// `chromInt()`.
    pub fn chrom_int(&self) -> Option<&ChromInterval> {
        self.chrom_int.as_ref()
    }
    /// `excludesamples()`.
    pub fn excludesamples(&self) -> Option<&Path> {
        self.excludesamples.as_deref()
    }
    /// `excludemarkers()`.
    pub fn excludemarkers(&self) -> Option<&Path> {
        self.excludemarkers.as_deref()
    }

    // phasing parameters
    /// `burnin()`.
    pub fn burnin(&self) -> i32 {
        self.burnin
    }
    /// `iterations()`.
    pub fn iterations(&self) -> i32 {
        self.iterations
    }
    /// `initial_lr()`.
    pub fn initial_lr(&self) -> f32 {
        self.initial_lr
    }
    /// `phase_states()`.
    pub fn phase_states(&self) -> i32 {
        self.phase_states
    }
    /// `step_scale()`.
    pub fn step_scale(&self) -> f32 {
        self.step_scale
    }
    /// `rare()`.
    pub fn rare(&self) -> f32 {
        self.rare
    }

    // imputation parameters
    /// `impute()`.
    pub fn impute(&self) -> bool {
        self.impute
    }
    /// `imp_states()`.
    pub fn imp_states(&self) -> i32 {
        self.imp_states
    }
    /// `imp_segment()`.
    pub fn imp_segment(&self) -> f32 {
        self.imp_segment
    }
    /// `imp_step()`.
    pub fn imp_step(&self) -> f32 {
        self.imp_step
    }
    /// `imp_nsteps()`.
    pub fn imp_nsteps(&self) -> i32 {
        self.imp_nsteps
    }
    /// `cluster()`.
    pub fn cluster(&self) -> f32 {
        self.cluster
    }
    /// `ap()`.
    pub fn ap(&self) -> bool {
        self.ap
    }
    /// `gp()`.
    pub fn gp(&self) -> bool {
        self.gp
    }

    // general parameters
    /// `em()`.
    pub fn em(&self) -> bool {
        self.em
    }
    /// `ne()`.
    pub fn ne(&self) -> f32 {
        self.ne
    }
    /// `err(int nHaps)` — the configured value, else the Li–Stephens default.
    pub fn err(&self, n_haps: i32) -> f32 {
        assert!(n_haps > 0, "{n_haps}");
        if self.err >= 0.0 {
            self.err
        } else {
            li_stephens_p_mismatch(n_haps)
        }
    }
    /// `window()`.
    pub fn window(&self) -> f32 {
        self.window
    }
    /// `window_markers()`.
    pub fn window_markers(&self) -> i32 {
        self.window_markers
    }
    /// `overlap()`.
    pub fn overlap(&self) -> f32 {
        self.overlap
    }
    /// `buffer()`.
    pub fn buffer(&self) -> f32 {
        self.buffer
    }
    /// `seed()`.
    pub fn seed(&self) -> i64 {
        self.seed
    }
    /// `nthreads()`.
    pub fn nthreads(&self) -> i32 {
        self.nthreads
    }
    /// `truth()`.
    pub fn truth(&self) -> Option<&Path> {
        self.truth.as_deref()
    }
}

/// `Par.liStephensPMismatch(int nHaps)`. Parity note (num-2): uses `ln`; Java `Math.log`
/// vs Rust `f64::ln` may differ by ~1 ULP (see docs/known-quirks.md).
pub fn li_stephens_p_mismatch(n_haps: i32) -> f32 {
    let theta = 1.0 / ((n_haps as f64).ln() + 0.5);
    (theta / (2.0 * (theta + n_haps as f64))) as f32
}

fn available_processors() -> i32 {
    std::thread::available_parallelism()
        .map(|n| n.get() as i32)
        .unwrap_or(1)
}

/// `String.valueOf(float)` for the in-range decimal values used in `usage()`. Java's
/// `Float.toString` always keeps a fractional digit (`40.0`, `2.0`), whereas Rust's `{}` drops a
/// trailing `.0` for integer-valued floats; this re-adds it. (Only used for the usage defaults,
/// which lie in `[1e-3, 1e7)`, so the shortest-decimal output of `{}` matches Java's digits.)
fn java_float_str(f: f32) -> String {
    let s = format!("{f}");
    if s.contains(['.', 'e', 'E', 'N', 'n']) {
        s
    } else {
        format!("{s}.0")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn touch(name: &str) -> PathBuf {
        let path = std::env::temp_dir().join(name);
        std::fs::File::create(&path)
            .unwrap()
            .write_all(b"x")
            .unwrap();
        path
    }

    #[test]
    fn java_float_str_keeps_fractional_digit() {
        assert_eq!(java_float_str(40.0), "40.0"); // Java: String.valueOf(40.0f)
        assert_eq!(java_float_str(2.0), "2.0");
        assert_eq!(java_float_str(0.005), "0.005");
    }

    #[test]
    fn usage_matches_java_layout() {
        let u = Par::usage();
        // argument lines keep the leading two-space indent (regression: `\`-continuation ate it)
        assert!(u.contains("\n  gt=<VCF file with GT FORMAT field>                 (required)\n"));
        // float defaults are rendered Java-style (40.0 / 2.0, not 40 / 2)
        assert!(
            u.contains("(default=40.0)"),
            "window default not Java-formatted"
        );
        assert!(
            u.contains("(default=2.0)"),
            "overlap default not Java-formatted"
        );
        assert!(u.starts_with("Usage: java -jar beagle.27Feb25.75f.jar [arguments]\n\n"));
        assert!(u.ends_with("(default: machine dependent)\n\n"));
    }

    #[test]
    fn parses_required_and_defaults() {
        let gt = touch("beagle_rs_par_gt.vcf");
        let args = vec![
            format!("gt={}", gt.display()),
            "out=result".to_string(),
            "nthreads=1".to_string(),
            "seed=99999".to_string(),
        ];
        let par = Par::new(&args);
        assert_eq!(par.gt(), gt.as_path());
        assert_eq!(par.out(), "result");
        assert_eq!(par.nthreads(), 1);
        assert!(!par.no_n_threads());
        assert_eq!(par.seed(), 99999);
        assert_eq!(par.burnin(), D_BURNIN);
        assert_eq!(par.iterations(), D_ITERATIONS);
        assert_eq!(par.phase_states(), D_PHASE_STATES);
        assert!(par.impute());
        assert_eq!(par.imp_states(), D_IMP_STATES);
        assert!(par.ref_file().is_none());
        assert!(par.ped().is_none());
        // default err is negative -> err(nHaps) returns the Li-Stephens default
        assert!(par.err(1000) > 0.0);
    }

    #[test]
    fn overrides_parsed() {
        let gt = touch("beagle_rs_par_gt2.vcf");
        let args = vec![
            format!("gt={}", gt.display()),
            "out=o".to_string(),
            "window=20.0".to_string(),
            "ne=50000".to_string(),
            "impute=false".to_string(),
            "gp=true".to_string(),
        ];
        let par = Par::new(&args);
        assert!((par.window() - 20.0).abs() < 1e-6);
        assert!((par.ne() - 50000.0).abs() < 1e-3);
        assert!(!par.impute());
        assert!(par.gp());
    }

    #[test]
    fn li_stephens_matches_formula() {
        let n = 2000;
        let theta = 1.0 / ((n as f64).ln() + 0.5);
        let expected = (theta / (2.0 * (theta + n as f64))) as f32;
        assert_eq!(li_stephens_p_mismatch(n), expected);
    }

    #[test]
    #[should_panic(expected = "value=0.6 > 0.5")]
    fn rejects_out_of_range_rare() {
        let gt = touch("beagle_rs_par_gt3.vcf");
        let args = vec![
            format!("gt={}", gt.display()),
            "out=o".to_string(),
            "rare=0.6".to_string(), // > 0.5 max -> panics in Validate
        ];
        let _ = Par::new(&args);
    }
}
