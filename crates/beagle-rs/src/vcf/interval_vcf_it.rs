//! Port of `vcf/IntervalVcfIt.java` — a `SampleFileIt` that skips leading records outside
//! a chromosome interval, then yields the contiguous run of records inside it (stopping at
//! the first record past the interval, assuming position-sorted input).

use std::path::Path;

use crate::beagleutil::ChromInterval;
use crate::blbutil::{FileIt, SampleFileIt};

use super::{GTRec, Samples};

type Rec = Box<dyn GTRec>;

/// Port of `vcf/IntervalVcfIt.java`.
pub struct IntervalVcfIt<I: SampleFileIt<Item = Rec>> {
    it: I,
    interval: ChromInterval,
    next: Option<Rec>,
}

fn read_first_record<I: SampleFileIt<Item = Rec>>(
    it: &mut I,
    interval: &ChromInterval,
) -> Option<Rec> {
    it.by_ref()
        .find(|candidate| interval.contains(candidate.marker()))
}

fn read_next_record<I: SampleFileIt<Item = Rec>>(
    it: &mut I,
    interval: &ChromInterval,
) -> Option<Rec> {
    match it.next() {
        Some(candidate) if interval.contains(candidate.marker()) => Some(candidate),
        _ => None,
    }
}

impl<I: SampleFileIt<Item = Rec>> IntervalVcfIt<I> {
    /// `new IntervalVcfIt(SampleFileIt it, ChromInterval chromInt)`.
    pub fn new(mut it: I, chrom_int: ChromInterval) -> Self {
        let first_record = read_first_record(&mut it, &chrom_int);
        let first_record = first_record.unwrap_or_else(|| {
            panic!(
                "No VCF records found in the specified interval.\n\
                 Check chromosome identifier and interval: {chrom_int}"
            )
        });
        IntervalVcfIt {
            it,
            interval: chrom_int,
            next: Some(first_record),
        }
    }
}

impl<I: SampleFileIt<Item = Rec>> Iterator for IntervalVcfIt<I> {
    type Item = Rec;

    fn next(&mut self) -> Option<Rec> {
        let current = self.next.take()?;
        self.next = read_next_record(&mut self.it, &self.interval);
        Some(current)
    }
}

impl<I: SampleFileIt<Item = Rec>> FileIt for IntervalVcfIt<I> {
    fn file(&self) -> Option<&Path> {
        self.it.file()
    }
}

impl<I: SampleFileIt<Item = Rec>> SampleFileIt for IntervalVcfIt<I> {
    fn samples(&self) -> &Samples {
        self.it.samples()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blbutil::InputIt;
    use crate::vcf::{to_basic_gt_rec, VcfIt};
    use std::io::{BufReader, Cursor};

    const VCF: &str = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0
chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1
chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t1|0
chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t0|0
chr1\t400\t.\tA\tT\t.\tPASS\t.\tGT\t1|1
";

    fn vcf_it() -> VcfIt<InputIt> {
        let it = InputIt::from_reader(Box::new(BufReader::new(Cursor::new(VCF.to_owned()))), None);
        VcfIt::create(it, to_basic_gt_rec)
    }

    #[test]
    fn yields_only_records_in_interval() {
        let interval = ChromInterval::new("chr1", 150, 350);
        let ivi = IntervalVcfIt::new(vcf_it(), interval);
        assert_eq!(ivi.samples().size(), 1);
        let positions: Vec<i32> = ivi.map(|r| r.marker().pos()).collect();
        assert_eq!(positions, vec![200, 300]);
    }

    #[test]
    #[should_panic(expected = "No VCF records found")]
    fn panics_when_interval_empty() {
        let interval = ChromInterval::new("chr2", 1, 10); // different chromosome
        let _ = IntervalVcfIt::new(vcf_it(), interval);
    }
}
