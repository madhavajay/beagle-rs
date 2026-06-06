//! Port of `vcf/FilterUtil.java` — static factories for marker and sample filters.

use crate::beagleutil::ChromInterval;
use crate::blbutil::{Filter, Utilities};
use std::collections::HashSet;
use std::path::Path;

use super::{marker_utils, Marker};

/// `markerFilter(File excludeMarkersFile)` — excludes markers whose id or `CHROM:POS`
/// matches a line of the file; accepts all if `None`.
pub fn marker_filter(exclude_markers_file: Option<&Path>) -> Filter<Marker> {
    match exclude_markers_file {
        None => Filter::accept_all(),
        Some(f) => exclude_id_filter(Utilities::id_set(Some(f))),
    }
}

/// `chromIntFilter(ChromInterval chromInterval)` — excludes markers outside the interval;
/// accepts all if `None`.
pub fn chrom_int_filter(chrom_interval: Option<ChromInterval>) -> Filter<Marker> {
    match chrom_interval {
        None => Filter::accept_all(),
        Some(ci) => Filter::predicate(move |marker: &Marker| ci.contains(marker)),
    }
}

/// `sampleFilter(File excludeSamplesFile)` — excludes samples listed in the file; accepts
/// all if `None`.
pub fn sample_filter(exclude_samples_file: Option<&Path>) -> Filter<String> {
    match exclude_samples_file {
        None => Filter::accept_all(),
        Some(f) => Filter::exclude(Utilities::id_set(Some(f))),
    }
}

/// `sampleFilter(File sampleFile, boolean includeFilter)` — include- or exclude-filter
/// over the ids in `sampleFile`; accepts all if `None`.
pub fn sample_filter_inc_exc(sample_file: Option<&Path>, include_filter: bool) -> Filter<String> {
    match sample_file {
        None => Filter::accept_all(),
        Some(f) => {
            let id_set = Utilities::id_set(Some(f));
            if include_filter {
                Filter::include(id_set)
            } else {
                Filter::exclude(id_set)
            }
        }
    }
}

/// `markerIsInSet(Marker marker, Set<String> set)` — true if any marker id or `CHROM:POS`
/// is in `set`.
pub fn marker_is_in_set(marker: &Marker, set: &HashSet<String>) -> bool {
    for id in marker_utils::ids(marker) {
        if set.contains(&id) {
            return true;
        }
    }
    let pos_id = format!("{}:{}", marker.chrom(), marker.pos());
    set.contains(&pos_id)
}

/// `excludeIdFilter(Collection<String> exclude)` — accepts markers whose id/`CHROM:POS` is
/// not in `exclude` (accepts all if `exclude` is empty).
pub fn exclude_id_filter(exclude: HashSet<String>) -> Filter<Marker> {
    if exclude.is_empty() {
        Filter::predicate(|_marker: &Marker| true)
    } else {
        Filter::predicate(move |marker: &Marker| !marker_is_in_set(marker, &exclude))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::MarkerParser;
    use std::io::Write;

    fn marker(line: &str) -> Marker {
        Marker::instance(line, &MarkerParser::new(true, true, true, true))
    }

    fn write_lines(name: &str, contents: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(name);
        std::fs::File::create(&path)
            .unwrap()
            .write_all(contents.as_bytes())
            .unwrap();
        path
    }

    #[test]
    fn marker_is_in_set_matches_id_or_pos() {
        let m = marker("chr1\t100\trs1;rs2\tA\tC\t.\tPASS\t.\tGT\t0|0");
        let by_id: HashSet<String> = ["rs2".to_string()].into_iter().collect();
        assert!(marker_is_in_set(&m, &by_id));
        let by_pos: HashSet<String> = ["chr1:100".to_string()].into_iter().collect();
        assert!(marker_is_in_set(&m, &by_pos));
        let neither: HashSet<String> = ["rsX".to_string(), "chr1:200".to_string()]
            .into_iter()
            .collect();
        assert!(!marker_is_in_set(&m, &neither));
    }

    #[test]
    fn exclude_id_filter_predicate() {
        let exclude: HashSet<String> = ["chr1:100".to_string()].into_iter().collect();
        let f = exclude_id_filter(exclude);
        assert!(!f.accept(&marker("chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0")));
        assert!(f.accept(&marker("chr1\t200\t.\tA\tC\t.\tPASS\t.\tGT\t0|0")));

        // empty exclude -> accept all
        let all = exclude_id_filter(HashSet::new());
        assert!(all.accept(&marker("chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0")));
    }

    #[test]
    fn marker_filter_from_file() {
        let path = write_lines("beagle_rs_excl_markers.txt", "rs9\nchr1:100\n");
        let f = marker_filter(Some(&path));
        assert!(!f.accept(&marker("chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0"))); // pos match
        assert!(!f.accept(&marker("chr2\t5\trs9\tA\tC\t.\tPASS\t.\tGT\t0|0"))); // id match
        assert!(f.accept(&marker("chr2\t6\trs8\tA\tC\t.\tPASS\t.\tGT\t0|0")));
        // None -> accept all
        assert!(marker_filter(None).accept(&marker("chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0")));
    }

    #[test]
    fn sample_filters() {
        let path = write_lines("beagle_rs_excl_samples.txt", "S1\nS3\n");
        let excl = sample_filter(Some(&path));
        assert!(!excl.accept(&"S1".to_string()));
        assert!(excl.accept(&"S2".to_string()));

        let inc = sample_filter_inc_exc(Some(&path), true);
        assert!(inc.accept(&"S1".to_string()));
        assert!(!inc.accept(&"S2".to_string()));

        assert!(sample_filter(None).accept(&"anything".to_string()));
    }

    #[test]
    fn chrom_int_filter_predicate() {
        let ci = ChromInterval::new("chrCIF", 100, 200);
        let f = chrom_int_filter(Some(ci));
        assert!(f.accept(&marker("chrCIF\t150\t.\tA\tC\t.\tPASS\t.\tGT\t0|0")));
        assert!(!f.accept(&marker("chrCIF\t250\t.\tA\tC\t.\tPASS\t.\tGT\t0|0")));
        assert!(chrom_int_filter(None).accept(&marker("chrZ\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0|0")));
    }
}
