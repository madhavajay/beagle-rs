//! Port of `phase/SampleSeg.java` — an immutable `[start, inclEnd]` marker segment for a
//! sample.

use std::cmp::Ordering;
use std::fmt;

use crate::beagleutil::IntInterval;
use crate::blbutil::consts;

/// Port of `phase/SampleSeg.java`.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct SampleSeg {
    sample: i32,
    start: i32,
    incl_end: i32,
}

impl SampleSeg {
    /// `new SampleSeg(int sample, int start, int inclEnd)`.
    pub fn new(sample: i32, start: i32, incl_end: i32) -> Self {
        assert!(start <= incl_end, "{start}");
        SampleSeg {
            sample,
            start,
            incl_end,
        }
    }

    /// `sample()`.
    pub fn sample(&self) -> i32 {
        self.sample
    }

    /// `SampleSeg.sampleComp()` — order by sample, then start, then inclEnd.
    pub fn sample_cmp(t1: &SampleSeg, t2: &SampleSeg) -> Ordering {
        t1.sample
            .cmp(&t2.sample)
            .then(t1.start.cmp(&t2.start))
            .then(t1.incl_end.cmp(&t2.incl_end))
    }
}

impl IntInterval for SampleSeg {
    fn start(&self) -> i32 {
        self.start
    }
    fn incl_end(&self) -> i32 {
        self.incl_end
    }
}

impl Ord for SampleSeg {
    /// `compareTo` — order by start, then inclEnd, then sample.
    fn cmp(&self, other: &Self) -> Ordering {
        self.start
            .cmp(&other.start)
            .then(self.incl_end.cmp(&other.incl_end))
            .then(self.sample.cmp(&other.sample))
    }
}

impl PartialOrd for SampleSeg {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for SampleSeg {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "[{}: {}{}{}]",
            self.sample,
            self.start,
            consts::HYPHEN,
            self.incl_end
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accessors_and_display() {
        let ss = SampleSeg::new(3, 10, 20);
        assert_eq!(ss.sample(), 3);
        assert_eq!(ss.start(), 10);
        assert_eq!(ss.incl_end(), 20);
        assert_eq!(ss.to_string(), "[3: 10-20]");
    }

    #[test]
    #[should_panic(expected = "10")]
    fn rejects_start_after_end() {
        let _ = SampleSeg::new(0, 10, 5);
    }

    #[test]
    fn compare_to_orders_by_start_end_sample() {
        let a = SampleSeg::new(5, 10, 20);
        let b = SampleSeg::new(0, 10, 20); // same start/end, smaller sample
        let c = SampleSeg::new(9, 10, 25); // same start, larger end
        let d = SampleSeg::new(0, 11, 12); // larger start
        assert!(b < a); // sample 0 < 5 (start/end equal)
        assert!(a < c); // end 20 < 25
        assert!(c < d); // start 10 < 11
    }

    #[test]
    fn sample_comp_orders_by_sample_first() {
        let a = SampleSeg::new(0, 100, 200);
        let b = SampleSeg::new(1, 1, 2);
        // sample_cmp: sample 0 < 1 even though a's start is larger
        assert_eq!(SampleSeg::sample_cmp(&a, &b), Ordering::Less);
        // natural Ord: a.start 100 > b.start 1 -> a is greater
        assert!(a > b);
    }
}
