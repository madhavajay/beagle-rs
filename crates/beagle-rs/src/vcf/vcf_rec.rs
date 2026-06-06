//! Port of `vcf/VcfRec.java` — a full VCF record `GTRec`: it parses genotypes via
//! `BasicGTRec`/`VcfRecGTParser` and provides delimiter-based access to QUAL/FILTER/INFO
//! and per-sample FORMAT subfields. The shared header is held via `Arc`.

use crate::blbutil::{consts, StringUtil};
use crate::ints::IntArray;
use std::collections::HashMap;
use std::sync::Arc;

use super::{BasicGTRec, GTRec, Marker, MarkerParser, Samples, VcfHeader, VcfRecGTParser};

const SAMPLE_OFFSET: i32 = 9;

/// Port of `vcf/VcfRec.java`.
pub struct VcfRec {
    vcf_header: Arc<VcfHeader>,
    vcf_record: String,
    delimiters: Vec<i32>,
    marker: Marker,
    format_fields: Vec<String>,
    format_map: HashMap<String, i32>,
    gt_rec: BasicGTRec,
}

/// `VcfRec.gtIndex(a1, a2)` — VCF genotype index for an allele pair.
pub fn gt_index(a1: i32, a2: i32) -> i32 {
    assert!(a1 >= 0, "a1<0: {a1}");
    assert!(a2 >= 0, "a2<0: {a2}");
    if a1 < a2 {
        (a2 * (a2 + 1)) / 2 + a1
    } else {
        (a1 * (a1 + 1)) / 2 + a2
    }
}

fn compute_delimiters(vcf_header: &VcfHeader, vcf_record: &str) -> Vec<i32> {
    let n_fields = vcf_header.n_header_fields();
    let bytes = vcf_record.as_bytes();
    let mut delimiters = vec![0i32; (n_fields + 1) as usize];
    delimiters[0] = -1;
    for j in 1..n_fields as usize {
        delimiters[j] = index_of_tab(bytes, delimiters[j - 1] + 1);
        if delimiters[j] == -1 {
            field_count_error(vcf_header, vcf_record);
        }
    }
    if index_of_tab(bytes, delimiters[n_fields as usize - 1] + 1) != -1 {
        field_count_error(vcf_header, vcf_record);
    }
    delimiters[n_fields as usize] = vcf_record.len() as i32;
    delimiters
}

fn index_of_tab(bytes: &[u8], from: i32) -> i32 {
    bytes[from as usize..]
        .iter()
        .position(|&b| b == b'\t')
        .map_or(-1, |o| from + o as i32)
}

fn field_count_error(vcf_header: &VcfHeader, vcf_record: &str) -> ! {
    let fields = StringUtil::get_fields(vcf_record, b'\t');
    panic!(
        "VCF header line has {} fields, but data line has {} fields\nFile source: {}",
        vcf_header.n_header_fields(),
        fields.len(),
        vcf_header.src()
    );
}

impl VcfRec {
    /// `new VcfRec(VcfHeader, String vcfRecord, MarkerParser)`.
    pub fn new(
        vcf_header: Arc<VcfHeader>,
        vcf_record: String,
        marker_parser: &MarkerParser,
    ) -> Self {
        let delimiters = compute_delimiters(&vcf_header, &vcf_record);
        let marker = Marker::instance(&vcf_record, marker_parser);
        let format = format_field(&vcf_record, &delimiters);
        let format_fields = parse_formats(format, &vcf_record);
        let format_map = format_to_index_map(&vcf_header, &vcf_record, &format_fields);
        assert!(
            format_map.contains_key("GT"),
            "Missing FORMAT/GT field: {vcf_record}"
        );
        let gt_rec = {
            let parser = VcfRecGTParser::new(&vcf_header, &vcf_record, marker_parser);
            BasicGTRec::from_parser(&parser)
        };
        VcfRec {
            vcf_header,
            vcf_record,
            delimiters,
            marker,
            format_fields,
            format_map,
            gt_rec,
        }
    }

    fn substr(&self, start: i32, end: i32) -> &str {
        &self.vcf_record[start as usize..end as usize]
    }

    fn format_subfield_end(&self, mut start: i32) -> i32 {
        let bytes = self.vcf_record.as_bytes();
        while (start as usize) < bytes.len() {
            let c = bytes[start as usize];
            if c == b':' || c == b'\t' {
                return start;
            }
            start += 1;
        }
        start
    }

    /// `qual()`.
    pub fn qual(&self) -> &str {
        self.substr(self.delimiters[5] + 1, self.delimiters[6])
    }

    /// `filter()`.
    pub fn filter(&self) -> &str {
        self.substr(self.delimiters[6] + 1, self.delimiters[7])
    }

    /// `info()`.
    pub fn info(&self) -> &str {
        self.substr(self.delimiters[7] + 1, self.delimiters[8])
    }

    /// `format()` — the FORMAT field, or `""` if absent.
    pub fn format(&self) -> &str {
        if self.delimiters.len() > 9 {
            self.substr(self.delimiters[8] + 1, self.delimiters[9])
        } else {
            ""
        }
    }

    /// `nFormatSubfields()`.
    pub fn n_format_subfields(&self) -> i32 {
        self.format_fields.len() as i32
    }

    /// `formatSubfield(int)`.
    pub fn format_subfield(&self, subfield_index: i32) -> &str {
        &self.format_fields[subfield_index as usize]
    }

    /// `hasFormat(String)`.
    pub fn has_format(&self, format_code: &str) -> bool {
        self.format_map.contains_key(format_code)
    }

    /// `formatIndex(String)`.
    pub fn format_index(&self, format_code: &str) -> i32 {
        self.format_map.get(format_code).copied().unwrap_or(-1)
    }

    /// `sampleData(int sample)` — the whole sample column.
    pub fn sample_data(&self, sample: i32) -> &str {
        let index = self.vcf_header.unfiltered_sample_index(sample);
        self.substr(
            self.delimiters[(index + SAMPLE_OFFSET) as usize] + 1,
            self.delimiters[(index + SAMPLE_OFFSET + 1) as usize],
        )
    }

    /// `sampleData(int sample, int subfieldIndex)`.
    pub fn sample_data_subfield(&self, sample: i32, subfield_index: i32) -> &str {
        assert!(
            subfield_index >= 0 && subfield_index < self.format_fields.len() as i32,
            "{subfield_index}"
        );
        let index = SAMPLE_OFFSET + self.vcf_header.unfiltered_sample_index(sample);
        let mut start = self.delimiters[index as usize] + 1;
        let len = self.vcf_record.len() as i32;
        let bytes = self.vcf_record.as_bytes();
        for _ in 0..subfield_index {
            let end = self.format_subfield_end(start);
            if end == len || bytes[end as usize] == b'\t' {
                return consts::MISSING_DATA_STRING;
            }
            start = end + 1;
        }
        let end = self.format_subfield_end(start);
        if end == start {
            consts::MISSING_DATA_STRING
        } else {
            self.substr(start, end)
        }
    }

    /// `sampleData(int sample, String formatCode)`.
    pub fn sample_data_code(&self, sample: i32, format_code: &str) -> &str {
        let format_index = *self
            .format_map
            .get(format_code)
            .unwrap_or_else(|| panic!("missing format data: {format_code}"));
        self.sample_data_subfield(sample, format_index)
    }

    /// `formatData(String formatCode)`.
    pub fn format_data(&self, format_code: &str) -> Vec<String> {
        let format_index = *self
            .format_map
            .get(format_code)
            .unwrap_or_else(|| panic!("missing format data: {format_code}"));
        (0..self.vcf_header.n_samples())
            .map(|j| self.sample_data_subfield(j, format_index).to_string())
            .collect()
    }

    /// `vcfHeader()`.
    pub fn vcf_header(&self) -> &VcfHeader {
        &self.vcf_header
    }
}

fn format_field<'a>(vcf_record: &'a str, delimiters: &[i32]) -> &'a str {
    if delimiters.len() > 9 {
        &vcf_record[(delimiters[8] + 1) as usize..delimiters[9] as usize]
    } else {
        ""
    }
}

fn parse_formats(formats: &str, vcf_record: &str) -> Vec<String> {
    assert!(
        formats != consts::MISSING_DATA_STRING && !formats.is_empty(),
        "missing format field: {vcf_record}"
    );
    let fields = StringUtil::get_fields(formats, b':');
    for f in &fields {
        assert!(
            !f.is_empty(),
            "missing format in format subfield list: {vcf_record}"
        );
    }
    fields.into_iter().map(|s| s.to_string()).collect()
}

fn format_to_index_map(
    vcf_header: &VcfHeader,
    vcf_record: &str,
    format_fields: &[String],
) -> HashMap<String, i32> {
    if vcf_header.n_samples() == 0 {
        return HashMap::new();
    }
    let mut map = HashMap::with_capacity(format_fields.len());
    for (j, f) in format_fields.iter().enumerate() {
        map.insert(f.clone(), j as i32);
    }
    if let Some(&gt) = map.get("GT") {
        assert!(gt == 0, "GT format is not first format: {vcf_record}");
    }
    map
}

impl IntArray for VcfRec {
    fn size(&self) -> i32 {
        2 * self.vcf_header.n_samples()
    }
    fn get(&self, hap: i32) -> i32 {
        self.gt_rec.get(hap)
    }
}

impl GTRec for VcfRec {
    fn samples(&self) -> &Samples {
        self.vcf_header.samples()
    }
    fn marker(&self) -> &Marker {
        &self.marker
    }
    fn is_phased_sample(&self, sample: i32) -> bool {
        self.gt_rec.is_phased_sample(sample)
    }
    fn is_phased(&self) -> bool {
        self.gt_rec.is_phased()
    }
}

impl std::fmt::Display for VcfRec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.vcf_record)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::HEADER_PREFIX;

    fn header() -> Arc<VcfHeader> {
        let lines = vec![
            "##fileformat=VCFv4.2".to_string(),
            format!("{HEADER_PREFIX}\tS0\tS1"),
        ];
        // both diploid
        Arc::new(VcfHeader::new_accept_all("src", &lines, &[true, true]))
    }

    #[test]
    fn parses_record_fields_and_genotypes() {
        let rec = "chr1\t100\trs1\tA\tC\t30\tPASS\tDP=9\tGT:DP\t0|1:8\t1/1:7".to_string();
        let r = VcfRec::new(header(), rec, &MarkerParser::new(true, true, true, true));
        assert_eq!(r.qual(), "30");
        assert_eq!(r.filter(), "PASS");
        assert_eq!(r.info(), "DP=9");
        assert_eq!(r.format(), "GT:DP");
        assert!(r.has_format("GT"));
        assert!(r.has_format("DP"));
        assert!(!r.has_format("GQ"));
        assert_eq!(r.sample_data(0), "0|1:8");
        assert_eq!(r.sample_data_code(0, "DP"), "8");
        assert_eq!(r.sample_data_code(1, "DP"), "7");
        assert_eq!(r.sample_data_subfield(0, 0), "0|1");
        assert_eq!(r.format_data("DP"), vec!["8", "7"]);
        // genotypes
        assert_eq!(r.size(), 4);
        assert_eq!(
            (0..4).map(|h| r.get(h)).collect::<Vec<_>>(),
            vec![0, 1, 1, 1]
        );
        assert!(r.is_phased_sample(0));
        assert!(!r.is_phased_sample(1));
        assert_eq!(r.marker().pos(), 100);
    }

    #[test]
    fn gt_index_values() {
        assert_eq!(gt_index(0, 0), 0);
        assert_eq!(gt_index(0, 1), 1);
        assert_eq!(gt_index(1, 1), 2);
        assert_eq!(gt_index(0, 2), 3);
        assert_eq!(gt_index(2, 0), 3); // symmetric
    }
}
