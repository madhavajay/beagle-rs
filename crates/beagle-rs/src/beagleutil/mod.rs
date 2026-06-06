//! Port of the Java `beagleutil` package: chromosome/sample id registries, the
//! interval interface, and PBWT update primitives.
//!
//! Ported so far: `ThreadSafeIndexer`, `ChromIds`, `SampleIds`, `IntInterval`.
//! (`ChromInterval`, `CompHapSegment`, `PbwtUpdater`, `PbwtDivUpdater` land next.)

mod chrom_ids;
mod chrom_interval;
mod comp_hap_segment;
mod pbwt_div_updater;
mod pbwt_updater;
mod sample_ids;
mod thread_safe_indexer;

pub use chrom_ids::ChromIds;
pub use chrom_interval::ChromInterval;
pub use comp_hap_segment::CompHapSegment;
pub use pbwt_div_updater::PbwtDivUpdater;
pub use pbwt_updater::PbwtUpdater;
pub use sample_ids::SampleIds;
pub use thread_safe_indexer::ThreadSafeIndexer;

use std::cmp::Ordering;

/// Port of `beagleutil/IntInterval.java` — an immutable interval of consecutive ints.
pub trait IntInterval {
    /// Start of the interval (inclusive).
    fn start(&self) -> i32;
    /// End of the interval (inclusive).
    fn incl_end(&self) -> i32;
}

/// `IntInterval.incEndComp()` — by increasing `start`, then increasing `inclEnd`.
pub fn cmp_inc_end(t1: &dyn IntInterval, t2: &dyn IntInterval) -> Ordering {
    if t1.start() != t2.start() {
        if t1.start() < t2.start() {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    } else if t1.incl_end() != t2.incl_end() {
        if t1.incl_end() < t2.incl_end() {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    } else {
        Ordering::Equal
    }
}

/// `IntInterval.decEndComp()` — by increasing `start`, then decreasing `inclEnd`.
pub fn cmp_dec_end(t1: &dyn IntInterval, t2: &dyn IntInterval) -> Ordering {
    if t1.start() != t2.start() {
        if t1.start() < t2.start() {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    } else if t1.incl_end() != t2.incl_end() {
        if t1.incl_end() > t2.incl_end() {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    } else {
        Ordering::Equal
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Iv(i32, i32);
    impl IntInterval for Iv {
        fn start(&self) -> i32 {
            self.0
        }
        fn incl_end(&self) -> i32 {
            self.1
        }
    }

    #[test]
    fn interval_comparators() {
        assert_eq!(cmp_inc_end(&Iv(1, 5), &Iv(2, 0)), Ordering::Less);
        assert_eq!(cmp_inc_end(&Iv(2, 3), &Iv(2, 9)), Ordering::Less);
        assert_eq!(cmp_inc_end(&Iv(2, 9), &Iv(2, 9)), Ordering::Equal);
        // dec end: same start, larger inclEnd sorts first
        assert_eq!(cmp_dec_end(&Iv(2, 9), &Iv(2, 3)), Ordering::Less);
        assert_eq!(cmp_dec_end(&Iv(1, 0), &Iv(2, 9)), Ordering::Less);
    }
}
