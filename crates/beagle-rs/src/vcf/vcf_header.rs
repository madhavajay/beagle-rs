//! Port of `vcf/VcfHeader.java` — the VCF meta-information lines + `#CHROM` header
//! line, with sample filtering. Its `to_string()` (used to write output) reproduces the
//! header byte-for-byte.

use crate::blbutil::{consts, Filter, StringUtil, Utilities};

use super::marker_utils;
use super::{Samples, VcfMetaInfo};

/// First nine tab-delimited fields of a VCF header line with sample data.
pub const HEADER_PREFIX: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
const SAMPLE_OFFSET: usize = 9;

/// Port of `vcf/VcfHeader.java`.
pub struct VcfHeader {
    src: String,
    meta_info_lines: Vec<VcfMetaInfo>,
    n_header_fields: i32,
    included_indices: Vec<i32>,
    samples: Samples,
}

impl VcfHeader {
    /// `new VcfHeader(src, lines, isDiploid)` — accepts all samples.
    pub fn new_accept_all(src: &str, lines: &[String], is_diploid: &[bool]) -> Self {
        Self::new(src, lines, is_diploid, &Filter::accept_all())
    }

    /// `new VcfHeader(src, lines, isDiploid, sampleFilter)`.
    pub fn new(
        src: &str,
        lines: &[String],
        is_diploid: &[bool],
        sample_filter: &Filter<String>,
    ) -> Self {
        check_header_lines(lines, src);
        let header_index = lines.len() - 1;
        let meta_info_lines: Vec<VcfMetaInfo> = lines[0..header_index]
            .iter()
            .map(|l| VcfMetaInfo::new(l))
            .collect();
        let header_fields = StringUtil::get_fields(&lines[header_index], b'\t');
        let n_header_fields = header_fields.len() as i32;
        let included_indices = included_indices(src, &header_fields, sample_filter);
        let samples = build_samples(&header_fields, is_diploid, &included_indices);
        VcfHeader {
            src: src.to_string(),
            meta_info_lines,
            n_header_fields,
            included_indices,
            samples,
        }
    }

    /// `VcfHeader.isDiploid(vcfRec)` — per-sample diploid flags from the first record.
    pub fn is_diploid(vcf_rec: &str) -> Vec<bool> {
        let bytes = vcf_rec.as_bytes();
        let start = (marker_utils::ninth_tab_pos(vcf_rec) + 1) as usize;
        let mut list = Vec::new();
        let mut no_allele_sep = true;
        for &c in &bytes[start..] {
            if c == b'\t' {
                list.push(!no_allele_sep);
                no_allele_sep = true;
            } else if c == b'/' || c == b'|' {
                no_allele_sep = false;
            }
        }
        list.push(!no_allele_sep);
        list
    }

    /// `src()`.
    pub fn src(&self) -> &str {
        &self.src
    }

    /// `nMetaInfoLines()`.
    pub fn n_meta_info_lines(&self) -> i32 {
        self.meta_info_lines.len() as i32
    }

    /// `metaInfoLine(int index)`.
    pub fn meta_info_line(&self, index: i32) -> &VcfMetaInfo {
        &self.meta_info_lines[index as usize]
    }

    /// `nHeaderFields()`.
    pub fn n_header_fields(&self) -> i32 {
        self.n_header_fields
    }

    /// `nUnfilteredSamples()`.
    pub fn n_unfiltered_samples(&self) -> i32 {
        (self.n_header_fields - SAMPLE_OFFSET as i32).max(0)
    }

    /// `unfilteredSampleIndex(int sample)`.
    pub fn unfiltered_sample_index(&self, sample: i32) -> i32 {
        self.included_indices[sample as usize]
    }

    /// `nSamples()`.
    pub fn n_samples(&self) -> i32 {
        self.samples.size()
    }

    /// `samples()`.
    pub fn samples(&self) -> &Samples {
        &self.samples
    }

    /// `sampleIds()`.
    pub fn sample_ids(&self) -> Vec<String> {
        self.samples.ids()
    }
}

fn check_header_lines(lines: &[String], src: &str) {
    assert!(
        !lines.is_empty(),
        "\n\nERROR: Missing the VCF meta information lines and the VCF header line\nVCF source: {src}\n"
    );
    let line = &lines[lines.len() - 1];
    assert!(
        line.starts_with(HEADER_PREFIX),
        "\n\nERROR: Missing the VCF header line.\nVCF source: {src}\n\
         The VCF header line must immediately follow the meta-information lines.\n\
         The fields of the VCF header line must be tab-delimited and begin with:\n{HEADER_PREFIX}\n"
    );
}

fn included_indices(src: &str, header_fields: &[&str], sample_filter: &Filter<String>) -> Vec<i32> {
    let n_unfiltered = header_fields.len().saturating_sub(SAMPLE_OFFSET);
    let mut included = Vec::with_capacity(n_unfiltered);
    for j in 0..n_unfiltered {
        if sample_filter.accept(&header_fields[SAMPLE_OFFSET + j].to_string()) {
            included.push(j as i32);
        }
    }
    if included.is_empty() {
        Utilities::exit(&format!(
            "\nError      :  All samples in the VCF file are excluded\nFile       :  {src}"
        ));
    }
    included
}

fn build_samples(header_fields: &[&str], is_diploid: &[bool], included_indices: &[i32]) -> Samples {
    let mut ids = Vec::with_capacity(included_indices.len());
    let mut restricted_is_diploid = Vec::with_capacity(included_indices.len());
    for &idx in included_indices {
        ids.push(header_fields[SAMPLE_OFFSET + idx as usize].to_string());
        restricted_is_diploid.push(is_diploid[idx as usize]);
    }
    Samples::new(&ids, &restricted_is_diploid)
}

impl std::fmt::Display for VcfHeader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for mi in &self.meta_info_lines {
            write!(f, "{mi}{}", consts::NL)?;
        }
        f.write_str(HEADER_PREFIX)?;
        for id in self.samples.ids() {
            write!(f, "{}{}", consts::TAB, id)?;
        }
        write!(f, "{}", consts::NL)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lines() -> Vec<String> {
        vec![
            "##fileformat=VCFv4.2".to_string(),
            "##source=test".to_string(),
            format!("{HEADER_PREFIX}\tS0\tS1\tS2"),
        ]
    }

    #[test]
    fn parses_header_and_samples() {
        let recs = "chr1\t1\t.\tA\tC\t.\t.\t.\tGT\t0|1\t0/1\t1"; // 3 samples: diploid,diploid,haploid
        let is_dip = VcfHeader::is_diploid(recs);
        assert_eq!(is_dip, vec![true, true, false]);
        let h = VcfHeader::new_accept_all("src", &lines(), &is_dip);
        assert_eq!(h.n_samples(), 3);
        assert_eq!(h.sample_ids(), vec!["S0", "S1", "S2"]);
        assert!(h.samples().is_diploid(0));
        assert!(!h.samples().is_diploid(2));
        assert_eq!(h.unfiltered_sample_index(2), 2);
        assert_eq!(h.n_meta_info_lines(), 2);
        // header round-trips byte-for-byte
        assert_eq!(
            h.to_string(),
            "##fileformat=VCFv4.2\n##source=test\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\tS1\tS2\n"
        );
    }

    #[test]
    fn sample_filter_excludes() {
        let is_dip = vec![true, true, false];
        let filter = Filter::exclude(["S1".to_string()]);
        let h = VcfHeader::new("src", &lines(), &is_dip, &filter);
        assert_eq!(h.n_samples(), 2);
        assert_eq!(h.sample_ids(), vec!["S0", "S2"]);
        // sample 1 (S2) maps to unfiltered index 2
        assert_eq!(h.unfiltered_sample_index(1), 2);
    }
}
