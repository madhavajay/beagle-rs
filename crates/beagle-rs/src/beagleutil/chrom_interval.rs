//! Port of `beagleutil/ChromInterval.java` — a chromosome interval `[start, end]` in
//! genome coordinates.

use super::{ChromIds, IntInterval};
use crate::blbutil::consts;
use crate::vcf::Marker;

/// Port of `beagleutil/ChromInterval.java`. Field order matches Java `compareTo`
/// (chromIndex, start, end) so the derived `Ord` is identical.
#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ChromInterval {
    chrom_index: i32,
    start: i32,
    end: i32,
}

fn is_valid_pos(s: &str, start_index: i32, end_index: i32) -> bool {
    if start_index == end_index {
        return true;
    }
    let bytes = s.as_bytes();
    let length = end_index - start_index;
    // beagleutil-1: Java checks `s.charAt(startIndex) == 0` (the NUL char, not '0'),
    // so this leading-zero guard is dead code. Preserved.
    if length > 1 && bytes[start_index as usize] == 0 {
        return false;
    }
    for j in start_index..end_index {
        if !bytes[j as usize].is_ascii_digit() {
            return false;
        }
    }
    true
}

impl ChromInterval {
    /// `new ChromInterval(String chrom, int start, int end)`.
    pub fn new(chrom: &str, start: i32, end: i32) -> Self {
        assert!(start <= end, "start={start} end={end}");
        ChromInterval {
            chrom_index: ChromIds::instance().get_index(chrom),
            start,
            end,
        }
    }

    /// `ChromInterval.parse(String)` — `[chrom]`, `[chrom]:`, `[chrom]:[s]-[e]`,
    /// `[chrom]:[s]-`, `[chrom]:-[e]`. Returns `None` for invalid input. Missing start
    /// → `i32::MIN`, missing end → `i32::MAX`.
    pub fn parse(str: &str) -> Option<ChromInterval> {
        let s = str.trim();
        let length = s.len() as i32;
        let mut start = i32::MIN;
        let mut end = i32::MAX;
        let chr_delim = s.rfind(consts::COLON).map_or(-1, |i| i as i32);
        let pos_delim = s.rfind(consts::HYPHEN).map_or(-1, |i| i as i32);
        if length == 0 {
            None
        } else if chr_delim == -1 {
            Some(ChromInterval::new(s, start, end))
        } else if chr_delim == length - 1 {
            Some(ChromInterval::new(&s[0..(length - 1) as usize], start, end))
        } else {
            if pos_delim == -1
                || pos_delim <= chr_delim
                || chr_delim == length - 2
                || !is_valid_pos(s, chr_delim + 1, pos_delim)
                || !is_valid_pos(s, pos_delim + 1, length)
            {
                return None;
            }
            if pos_delim > chr_delim + 1 {
                start = s[(chr_delim + 1) as usize..pos_delim as usize]
                    .parse()
                    .expect("valid start position");
            }
            if length > pos_delim + 1 {
                end = s[(pos_delim + 1) as usize..length as usize]
                    .parse()
                    .expect("valid end position");
            }
            if start > end {
                return None;
            }
            Some(ChromInterval::new(&s[0..chr_delim as usize], start, end))
        }
    }

    /// `new ChromInterval(Marker start, Marker end)`.
    pub fn from_markers(start: &Marker, end: &Marker) -> Self {
        assert!(
            start.chrom_index() == end.chrom_index(),
            "start.chromIndex() != end.chromIndex()"
        );
        assert!(
            start.pos() >= 0 && start.pos() <= end.pos(),
            "start={start} end={end}"
        );
        ChromInterval {
            chrom_index: start.chrom_index(),
            start: start.pos(),
            end: end.pos(),
        }
    }

    /// `contains(Marker marker)` — true iff the marker lies in this interval.
    pub fn contains(&self, marker: &Marker) -> bool {
        let pos = marker.pos();
        marker.chrom_index() == self.chrom_index && self.start <= pos && pos <= self.end
    }

    /// `chromIndex()`.
    pub fn chrom_index(&self) -> i32 {
        self.chrom_index
    }

    /// `chrom()`.
    pub fn chrom(&self) -> String {
        ChromIds::instance().id(self.chrom_index)
    }

    /// `ChromInterval.overlap(a, b)`.
    pub fn overlap(a: &ChromInterval, b: &ChromInterval) -> bool {
        a.chrom_index == b.chrom_index && a.start <= b.end && b.start <= a.end
    }

    /// `ChromInterval.merge(a, b)` — union of overlapping intervals.
    pub fn merge(a: &ChromInterval, b: &ChromInterval) -> ChromInterval {
        assert!(
            ChromInterval::overlap(a, b),
            "non-overlapping intervals: {a} {b}"
        );
        ChromInterval::new(&a.chrom(), a.start.min(b.start), a.end.max(b.end))
    }
}

impl IntInterval for ChromInterval {
    fn start(&self) -> i32 {
        self.start
    }
    fn incl_end(&self) -> i32 {
        self.end
    }
}

impl std::fmt::Display for ChromInterval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", ChromIds::instance().id(self.chrom_index))?;
        if self.start > i32::MIN || self.end < i32::MAX {
            write!(f, "{}", consts::COLON)?;
            if self.start > i32::MIN {
                write!(f, "{}", self.start)?;
            }
            write!(f, "{}", consts::HYPHEN)?;
            if self.end < i32::MAX {
                write!(f, "{}", self.end)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_forms() {
        let c = ChromInterval::parse("chrPARSE1").unwrap();
        assert_eq!(c.chrom(), "chrPARSE1");
        assert_eq!(c.start(), i32::MIN);
        assert_eq!(c.incl_end(), i32::MAX);

        let c = ChromInterval::parse("chrPARSE2:100-200").unwrap();
        assert_eq!((c.start(), c.incl_end()), (100, 200));
        let c = ChromInterval::parse("chrPARSE2:100-").unwrap();
        assert_eq!((c.start(), c.incl_end()), (100, i32::MAX));
        let c = ChromInterval::parse("chrPARSE2:-200").unwrap();
        assert_eq!((c.start(), c.incl_end()), (i32::MIN, 200));
        let c = ChromInterval::parse("chrPARSE2:").unwrap();
        assert_eq!((c.start(), c.incl_end()), (i32::MIN, i32::MAX));

        assert!(ChromInterval::parse("").is_none());
        assert!(ChromInterval::parse("chrPARSE2:200-100").is_none());
        assert!(ChromInterval::parse("chrPARSE2:a-200").is_none());
    }

    #[test]
    fn display_forms() {
        assert_eq!(
            ChromInterval::new("chrDISP", 100, 200).to_string(),
            "chrDISP:100-200"
        );
        assert_eq!(
            ChromInterval::new("chrDISP", i32::MIN, i32::MAX).to_string(),
            "chrDISP"
        );
        assert_eq!(
            ChromInterval::new("chrDISP", 100, i32::MAX).to_string(),
            "chrDISP:100-"
        );
        assert_eq!(
            ChromInterval::new("chrDISP", i32::MIN, 200).to_string(),
            "chrDISP:-200"
        );
    }

    #[test]
    fn overlap_and_merge() {
        let a = ChromInterval::new("chrOV", 100, 200);
        let b = ChromInterval::new("chrOV", 150, 300);
        let c = ChromInterval::new("chrOV", 400, 500);
        assert!(ChromInterval::overlap(&a, &b));
        assert!(!ChromInterval::overlap(&a, &c));
        let m = ChromInterval::merge(&a, &b);
        assert_eq!((m.start(), m.incl_end()), (100, 300));
    }
}
