//! Port of `vcf/VcfWriter.java` — static methods for writing data in VCF 4.2 format.
//!
//! The Java class is a `final` holder of static methods with a private constructor; the
//! Rust port mirrors it as a set of free functions in this module.
//!
//! Parity note: `now()` (the `##filedate=` field) reads the wall clock and is therefore
//! the one documented non-reproducible field — see `docs/known-quirks.md`. Java uses
//! `Calendar.getInstance()` (the JVM default timezone); the Rust port computes the civil
//! date from `SystemTime` in UTC. This field is excluded from byte-for-byte comparisons.

use crate::blbutil::consts;
use std::io::{self, Write};
use std::time::{SystemTime, UNIX_EPOCH};

use super::{GTRec, Marker, Markers, Samples, VcfRecBuilder, GT};

const FILE_FORMAT: &str = "##fileformat=VCFv4.2";

const AF_INFO: &str = "##INFO=<ID=AF,Number=A,Type=Float,\
Description=\"Estimated ALT Allele Frequencies\">";
const DR2_INFO: &str = "##INFO=<ID=DR2,Number=A,Type=Float,\
Description=\"Dosage R-Squared: estimated squared correlation between \
estimated REF dose [P(RA) + 2*P(RR)] and true REF dose\">";
const IMP_INFO: &str = "##INFO=<ID=IMP,Number=0,Type=Flag,\
Description=\"Imputed marker\">";

const GT_FORMAT: &str = "##FORMAT=<ID=GT,Number=1,Type=String,\
Description=\"Genotype\">";
const DS_FORMAT: &str = "##FORMAT=<ID=DS,Number=A,Type=Float,\
Description=\"estimated ALT dose [P(RA) + 2*P(AA)]\">";
const AP1_FORMAT: &str = "##FORMAT=<ID=AP1,Number=A,Type=Float,\
Description=\"estimated ALT dose on first haplotype\">";
const AP2_FORMAT: &str = "##FORMAT=<ID=AP2,Number=A,Type=Float,\
Description=\"estimated ALT dose on second haplotype\">";
const GL_FORMAT: &str = "##FORMAT=<ID=GL,Number=G,Type=Float,\
Description=\"Log10-scaled Genotype Likelihood\">";
const GP_FORMAT: &str = "##FORMAT=<ID=GP,Number=G,Type=Float,\
Description=\"Estimated Genotype Probability\">";

// Mirrors the Java `SHORT_CHROM_PREFIX` field, which is itself only used to build
// `LONG_CHROM_PREFIX`. Kept for fidelity; `LONG_CHROM_PREFIX` is what gets written.
#[allow(dead_code)]
const SHORT_CHROM_PREFIX: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

const LONG_CHROM_PREFIX: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

/// `writeMetaLinesGT(String[] sampleIds, String source, PrintWriter out)` — only the GT
/// FORMAT subfield is described.
pub fn write_meta_lines_gt(
    sample_ids: &[String],
    source: Option<&str>,
    out: &mut dyn Write,
) -> io::Result<()> {
    write_meta_lines(sample_ids, source, false, false, false, false, out)
}

/// `writeMetaLines(String[] sampleIds, String source, boolean ds, boolean ap, boolean gp,
/// boolean gl, PrintWriter out)`.
#[allow(clippy::too_many_arguments)]
pub fn write_meta_lines(
    sample_ids: &[String],
    source: Option<&str>,
    ds: bool,
    ap: bool,
    gp: bool,
    gl: bool,
    out: &mut dyn Write,
) -> io::Result<()> {
    out.write_all(FILE_FORMAT.as_bytes())?;
    out.write_all(consts::NL.as_bytes())?;
    out.write_all(b"##filedate=")?;
    out.write_all(now().as_bytes())?;
    out.write_all(consts::NL.as_bytes())?;
    if let Some(src) = source {
        out.write_all(b"##source=\"")?;
        out.write_all(src.as_bytes())?;
        out.write_all(b"\"")?;
        out.write_all(consts::NL.as_bytes())?;
    }
    if ds {
        writeln_str(out, AF_INFO)?;
        writeln_str(out, DR2_INFO)?;
        writeln_str(out, IMP_INFO)?;
    }
    writeln_str(out, GT_FORMAT)?;
    if ds {
        writeln_str(out, DS_FORMAT)?;
    }
    if ap {
        writeln_str(out, AP1_FORMAT)?;
        writeln_str(out, AP2_FORMAT)?;
    }
    if gp {
        writeln_str(out, GP_FORMAT)?;
    }
    if gl {
        writeln_str(out, GL_FORMAT)?;
    }
    out.write_all(LONG_CHROM_PREFIX.as_bytes())?;
    for id in sample_ids {
        out.write_all(&[consts::TAB as u8])?;
        out.write_all(id.as_bytes())?;
    }
    out.write_all(consts::NL.as_bytes())
}

fn writeln_str(out: &mut dyn Write, s: &str) -> io::Result<()> {
    out.write_all(s.as_bytes())?;
    out.write_all(consts::NL.as_bytes())
}

/// `now()` — current date as `yyyyMMdd`. See the parity note at the top of this module:
/// this is the wall-clock exception and is computed in UTC.
fn now() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let days = (secs / 86400) as i64;
    let (y, m, d) = civil_from_days(days);
    format!("{y:04}{m:02}{d:02}")
}

/// Howard Hinnant's `civil_from_days`: days since 1970-01-01 → (year, month, day).
fn civil_from_days(days: i64) -> (i64, i64, i64) {
    let z = days + 719468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = z - era * 146097; // [0, 146096]
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365; // [0, 399]
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100); // [0, 365]
    let mp = (5 * doy + 2) / 153; // [0, 11]
    let d = doy - (153 * mp + 2) / 5 + 1; // [1, 31]
    let m = if mp < 10 { mp + 3 } else { mp - 9 }; // [1, 12]
    (if m <= 2 { y + 1 } else { y }, m, d)
}

/// `appendRecords(GT phasedTarg, int start, int end, PrintWriter out)` — phased target
/// genotypes for markers `[start, end)`.
pub fn append_records_gt(
    phased_targ: &dyn GT,
    start: i32,
    end: i32,
    out: &mut dyn Write,
) -> io::Result<()> {
    assert!(start <= end, "start={start} end={end}");
    assert!(phased_targ.is_phased(), "unphased genotypes");
    let samples = phased_targ.samples();
    let markers = phased_targ.markers();
    let mut rec_builders = rec_builders_range(markers, samples.size(), start, end);
    let n = phased_targ.n_samples();
    for s in 0..n {
        if samples.is_diploid(s) {
            let h1 = s << 1;
            let h2 = h1 | 0b1;
            for m in start..end {
                let a1 = phased_targ.allele(m, h1);
                let a2 = phased_targ.allele(m, h2);
                rec_builders[(m - start) as usize].add_sample_data_diploid(a1, a2);
            }
        } else {
            for m in start..end {
                let a1 = phased_targ.allele(m, s << 1);
                rec_builders[(m - start) as usize].add_sample_data_haploid(a1);
            }
        }
    }
    for vrb in &rec_builders {
        vrb.write_rec(out)?;
    }
    Ok(())
}

/// `appendRecords(Samples samples, GTRec[] recs, PrintWriter out)`.
pub fn append_records_recs(
    samples: &Samples,
    recs: &[&dyn GTRec],
    out: &mut dyn Write,
) -> io::Result<()> {
    let mut rec_builders = rec_builders_recs(samples, recs);
    let n_samples = samples.size();
    for s in 0..n_samples {
        if samples.is_diploid(s) {
            let h1 = s << 1;
            let h2 = h1 | 0b1;
            for (j, vrb) in rec_builders.iter_mut().enumerate() {
                let a1 = recs[j].get(h1);
                let a2 = recs[j].get(h2);
                vrb.add_sample_data_diploid(a1, a2);
            }
        } else {
            for (j, vrb) in rec_builders.iter_mut().enumerate() {
                let a1 = recs[j].get(s << 1);
                vrb.add_sample_data_haploid(a1);
            }
        }
    }
    for vrb in &rec_builders {
        vrb.write_rec(out)?;
    }
    Ok(())
}

fn rec_builders_recs(samples: &Samples, recs: &[&dyn GTRec]) -> Vec<VcfRecBuilder> {
    let n_samples = samples.size();
    recs.iter()
        .map(|rec| {
            assert!(rec.is_phased(), "unphased genotypes");
            assert!(rec.samples() == samples, "inconsistent samples");
            VcfRecBuilder::new(rec.marker().clone(), n_samples)
        })
        .collect()
}

fn rec_builders_range(
    markers: &Markers,
    n_samples: i32,
    start: i32,
    end: i32,
) -> Vec<VcfRecBuilder> {
    (start..end)
        .map(|m| VcfRecBuilder::new(markers.marker(m).clone(), n_samples))
        .collect()
}

/// `printFixedFieldsGT(Marker marker, PrintWriter out)` — the first 9 VCF record fields
/// with QUAL=`.`, FILTER=`PASS`, INFO=`.`, FORMAT=`GT` (the marker's own QUAL/FILTER/INFO
/// are intentionally overridden here).
pub fn print_fixed_fields_gt(marker: &Marker, out: &mut dyn Write) -> io::Result<()> {
    write!(out, "{marker}")?;
    out.write_all(&[consts::TAB as u8])?;
    out.write_all(&[consts::MISSING_DATA_CHAR as u8])?; // QUAL
    out.write_all(&[consts::TAB as u8])?;
    out.write_all(b"PASS")?; // FILTER
    out.write_all(&[consts::TAB as u8])?;
    out.write_all(&[consts::MISSING_DATA_CHAR as u8])?; // INFO
    out.write_all(&[consts::TAB as u8])?;
    out.write_all(b"GT")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ints::IntArray;
    use crate::vcf::MarkerParser;

    fn marker(line: &str) -> Marker {
        Marker::instance(line, &MarkerParser::new(true, true, true, true))
    }

    #[test]
    fn meta_lines_gt_structure() {
        let ids: Vec<String> = vec!["s0".into(), "s1".into()];
        let mut out = Vec::new();
        write_meta_lines_gt(&ids, None, &mut out).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[0], "##fileformat=VCFv4.2");
        assert!(lines[1].starts_with("##filedate="));
        let date = &lines[1]["##filedate=".len()..];
        assert_eq!(date.len(), 8);
        assert!(date.bytes().all(|b| b.is_ascii_digit()));
        assert_eq!(
            lines[2],
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        );
        assert_eq!(
            lines[3],
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1"
        );
        assert_eq!(lines.len(), 4);
        assert!(s.ends_with('\n'));
    }

    #[test]
    fn meta_lines_with_source_and_all_formats() {
        let ids: Vec<String> = vec!["a".into()];
        let mut out = Vec::new();
        write_meta_lines(&ids, Some("beagle"), true, true, true, true, &mut out).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[0], "##fileformat=VCFv4.2");
        assert!(lines[1].starts_with("##filedate="));
        assert_eq!(lines[2], "##source=\"beagle\"");
        // ds=true -> AF, DR2, IMP info lines
        assert!(lines[3].starts_with("##INFO=<ID=AF"));
        assert!(lines[4].starts_with("##INFO=<ID=DR2"));
        assert!(lines[5].starts_with("##INFO=<ID=IMP"));
        assert!(lines[6].starts_with("##FORMAT=<ID=GT"));
        assert!(lines[7].starts_with("##FORMAT=<ID=DS"));
        assert!(lines[8].starts_with("##FORMAT=<ID=AP1"));
        assert!(lines[9].starts_with("##FORMAT=<ID=AP2"));
        assert!(lines[10].starts_with("##FORMAT=<ID=GP"));
        assert!(lines[11].starts_with("##FORMAT=<ID=GL"));
        assert!(lines[12].starts_with("#CHROM"));
        assert!(lines[12].ends_with("\ta"));
    }

    #[test]
    fn print_fixed_fields_overrides_marker_fields() {
        // Marker has QUAL=50 FILTER=q10 INFO=foo, but printFixedFieldsGT overrides them.
        let m = marker("chr1\t100\trs1\tA\tC\t50\tq10\tfoo\tGT\t0|1");
        let mut out = Vec::new();
        print_fixed_fields_gt(&m, &mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "chr1\t100\trs1\tA\tC\t.\tPASS\t.\tGT"
        );
    }

    // Minimal phased GTRec for exercising append_records_recs.
    struct TestRec {
        marker: Marker,
        samples: Samples,
        alleles: Vec<i32>,
    }
    impl IntArray for TestRec {
        fn size(&self) -> i32 {
            self.alleles.len() as i32
        }
        fn get(&self, index: i32) -> i32 {
            self.alleles[index as usize]
        }
    }
    impl GTRec for TestRec {
        fn samples(&self) -> &Samples {
            &self.samples
        }
        fn marker(&self) -> &Marker {
            &self.marker
        }
        fn is_phased_sample(&self, _sample: i32) -> bool {
            true
        }
        fn is_phased(&self) -> bool {
            true
        }
    }

    #[test]
    fn append_records_writes_record_lines() {
        let samples = Samples::new(&["s0".into(), "s1".into()], &[true, true]);
        let r0 = TestRec {
            marker: marker("chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1"),
            samples: samples.clone(),
            alleles: vec![0, 1, 1, 0], // s0=0|1, s1=1|0
        };
        let r1 = TestRec {
            marker: marker("chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|0"),
            samples: samples.clone(),
            alleles: vec![1, 1, 0, 1], // s0=1|1, s1=0|1
        };
        let recs: Vec<&dyn GTRec> = vec![&r0, &r1];
        let mut out = Vec::new();
        append_records_recs(&samples, &recs, &mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0\n\
             chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t1|1\t0|1\n"
        );
    }

    #[test]
    fn civil_from_days_known_dates() {
        assert_eq!(civil_from_days(0), (1970, 1, 1));
        // 2000-01-01 is 10957 days after the epoch.
        assert_eq!(civil_from_days(10957), (2000, 1, 1));
        // 2021-03-01.
        assert_eq!(civil_from_days(18687), (2021, 3, 1));
    }
}
