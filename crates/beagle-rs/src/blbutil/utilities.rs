//! Port of `blbutil/Utilities.java` (the package-independent parts). The remaining
//! IO-dependent helpers (`duoPrint*`) and runtime/time helpers (`timeStamp`,
//! `printMemoryUse`, `commandLine`) are ported alongside `main`.

use crate::blbutil::{InputIt, StringUtil};
use crate::jdk::Random;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::path::Path;

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
