//! Port of `blbutil/Utilities.java` (the package-independent parts) plus the IO/runtime
//! helpers used by `main` (`duoPrint*`, `timeStamp`, `commandLine`).
//!
//! Parity note: `timeStamp()` reads the wall clock and is a documented `.log` normalization
//! field (run date), like `vcf::VcfWriter`'s `filedate`. `commandLine()` echoes the invocation;
//! the Rust binary has no JVM `-Xmx` heap setting, so that `-Xmx<N>m` token is omitted (another
//! `.log` line that differs from the `java -jar` reference and is normalized in comparisons).

use crate::blbutil::{consts, InputIt, StringUtil};
use crate::jdk::Random;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::io::Write;
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

/// Port of `blbutil/Utilities.java` (static methods).
pub struct Utilities;

impl Utilities {
    /// `shuffle(int[] ia, Random random)` — Fisher–Yates using `nextInt(bound)`.
    pub fn shuffle(ia: &mut [i32], random: &mut Random) {
        let n = ia.len();
        for j in 0..n {
            let x = random.next_int_bound((n - j) as i32) as usize;
            ia.swap(j, j + x);
        }
    }

    /// `shuffle(int[] ia, int nElements, Random random)`.
    pub fn shuffle_n(ia: &mut [i32], n_elements: i32, random: &mut Random) {
        let n = ia.len();
        for j in 0..n_elements as usize {
            let x = random.next_int_bound((n - j) as i32) as usize;
            ia.swap(j, j + x);
        }
    }

    /// `elapsedNanos(long nanoseconds)` — "H hours M minutes S seconds".
    pub fn elapsed_nanos(nanoseconds: i64) -> String {
        let mut seconds = (nanoseconds as f64 / 1_000_000_000.0).round() as i64;
        let mut sb = String::with_capacity(80);
        if seconds >= 3600 {
            let hours = seconds / 3600;
            sb.push_str(&hours.to_string());
            sb.push_str(if hours == 1 { " hour " } else { " hours " });
            seconds %= 3600;
        }
        if seconds >= 60 {
            let minutes = seconds / 60;
            sb.push_str(&minutes.to_string());
            sb.push_str(if minutes == 1 {
                " minute "
            } else {
                " minutes "
            });
            seconds %= 60;
        }
        sb.push_str(&seconds.to_string());
        sb.push_str(if seconds == 1 { " second" } else { " seconds" });
        sb
    }

    /// `commandLine(String program, String[] args)` — the multi-line command echo for the log.
    /// The JVM `-Xmx<N>m` token has no Rust equivalent and is omitted (see the module note).
    pub fn command_line(program: &str, args: &[String]) -> String {
        let mut sb = String::with_capacity(args.len() * 20);
        sb.push_str(consts::NL);
        sb.push_str("Command line: java");
        sb.push_str(" -jar ");
        sb.push_str(program);
        sb.push_str(consts::NL);
        for arg in args {
            sb.push_str("  ");
            sb.push_str(arg);
            sb.push_str(consts::NL);
        }
        sb
    }

    /// Seconds since the Unix epoch for wall-clock output fields (VCF `##filedate`, `.log`
    /// start/end timestamps). Honors `SOURCE_DATE_EPOCH` (the reproducible-builds standard) when
    /// it is set to a valid non-negative integer, so output can be made deterministic; otherwise
    /// reads the system clock. This is a Rust-port convenience — Java always reads the clock —
    /// and only affects fields already excluded from byte-for-byte comparison.
    pub fn wall_clock_secs() -> u64 {
        if let Ok(s) = std::env::var("SOURCE_DATE_EPOCH") {
            if let Ok(epoch) = s.trim().parse::<u64>() {
                return epoch;
            }
        }
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_secs())
            .unwrap_or(0)
    }

    /// `timeStamp()` — current UTC time as `hh:mm a 'UTC on' dd MMM yyyy` (wall-clock; a
    /// documented `.log` normalization field).
    pub fn time_stamp() -> String {
        let secs = Self::wall_clock_secs();
        let days = (secs / 86400) as i64;
        let sod = (secs % 86400) as i64; // seconds of day (UTC)
        let hour24 = sod / 3600;
        let minute = (sod % 3600) / 60;
        let am_pm = if hour24 < 12 { "AM" } else { "PM" };
        let mut hour12 = hour24 % 12;
        if hour12 == 0 {
            hour12 = 12;
        }
        let (y, m, d) = civil_from_days(days);
        const MONTHS: [&str; 12] = [
            "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
        ];
        format!(
            "{hour12:02}:{minute:02} {am_pm} UTC on {d:02} {mon} {y:04}",
            mon = MONTHS[(m - 1) as usize]
        )
    }

    /// `duoPrint(PrintWriter out, CharSequence s)` — print `s` to stdout and to `out`.
    pub fn duo_print(out: &mut dyn Write, s: &str) {
        print!("{s}");
        let _ = out.write_all(s.as_bytes());
    }

    /// `duoPrintln(PrintWriter out, CharSequence s)` — `println` `s` to stdout and to `out`.
    pub fn duo_println(out: &mut dyn Write, s: &str) {
        println!("{s}");
        let _ = out.write_all(s.as_bytes());
        let _ = out.write_all(b"\n");
    }

    /// `arrayToMap(E[] array)` — element → index map; panics on a duplicate element.
    pub fn array_to_map<E: Eq + Hash + Clone>(array: &[E]) -> HashMap<E, i32> {
        let mut map = HashMap::with_capacity(array.len());
        for (j, e) in array.iter().enumerate() {
            if map.insert(e.clone(), j as i32).is_some() {
                panic!("duplicate element in array");
            }
        }
        map
    }

    /// `commonIndices(E[] a1, E[] a2)` — `[a1-indices, a2-indices]` of common elements,
    /// in `a1` order. Panics if either array has a duplicate element.
    pub fn common_indices<E: Eq + Hash + Clone>(a1: &[E], a2: &[E]) -> [Vec<i32>; 2] {
        let map1 = Self::array_to_map(a1);
        let map2 = Self::array_to_map(a2);
        let mut idx0 = Vec::new();
        let mut idx1 = Vec::new();
        for id in a1 {
            if let Some(&j2) = map2.get(id) {
                idx0.push(map1[id]);
                idx1.push(j2);
            }
        }
        [idx0, idx1]
    }

    /// `idSet(File file)` — the set of trimmed, non-empty single-field lines of `file`
    /// (`None` → empty set). Panics if the file is missing, a directory, or any line has
    /// more than one white-space-delimited field.
    pub fn id_set(file: Option<&Path>) -> HashSet<String> {
        let mut id_set = HashSet::new();
        let Some(file) = file else {
            return id_set;
        };
        assert!(file.exists(), "file does not exist: {}", file.display());
        assert!(!file.is_dir(), "file is a directory: {}", file.display());
        for line in InputIt::from_gzip_file(file) {
            let line = line.trim();
            if !line.is_empty() {
                assert!(
                    StringUtil::count_fields_ws(line) <= 1,
                    "Line has more than one white-space delimited field (file: '{}'; line: '{}')",
                    file.display(),
                    line
                );
                id_set.insert(line.to_string());
            }
        }
        id_set
    }

    /// `exit(String s)` — print to stderr and terminate with exit code 1.
    pub fn exit(s: &str) -> ! {
        eprintln!();
        eprintln!("{s}");
        eprintln!();
        eprintln!("Terminating program.");
        std::process::exit(1);
    }
}

/// Howard Hinnant's `civil_from_days`: days since 1970-01-01 → (year, month, day).
fn civil_from_days(days: i64) -> (i64, i64, i64) {
    let z = days + 719468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    (if m <= 2 { y + 1 } else { y }, m, d)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shuffle_is_a_permutation() {
        let mut ia: Vec<i32> = (0..50).collect();
        let mut r = Random::new(7);
        Utilities::shuffle(&mut ia, &mut r);
        let mut sorted = ia.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..50).collect::<Vec<_>>());
    }

    #[test]
    fn elapsed_nanos_formats() {
        assert_eq!(Utilities::elapsed_nanos(0), "0 seconds");
        assert_eq!(Utilities::elapsed_nanos(1_000_000_000), "1 second");
        assert_eq!(
            Utilities::elapsed_nanos(61_000_000_000),
            "1 minute 1 second"
        );
        assert_eq!(
            Utilities::elapsed_nanos(3_661_000_000_000),
            "1 hour 1 minute 1 second"
        );
        assert_eq!(
            Utilities::elapsed_nanos(7_322_000_000_000),
            "2 hours 2 minutes 2 seconds"
        );
    }

    #[test]
    fn common_indices_and_array_to_map() {
        let a1 = ["a", "b", "c", "d"];
        let a2 = ["x", "c", "a", "y"];
        let [i0, i1] = Utilities::common_indices(&a1, &a2);
        // common elements in a1 order: "a" (a1[0],a2[2]), "c" (a1[2],a2[1])
        assert_eq!(i0, vec![0, 2]);
        assert_eq!(i1, vec![2, 1]);
    }

    #[test]
    #[should_panic]
    fn array_to_map_rejects_duplicates() {
        let _ = Utilities::array_to_map(&["a", "a"]);
    }
}
