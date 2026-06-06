//! Port of `vcf/VcfIt.java` — an iterator whose items are `GTRec`s parsed from the data
//! lines of a VCF file.
//!
//! Java buffers records and parses them in parallel batches; the buffer size derives from
//! `Runtime.maxMemory()`. Buffering only affects granularity, never the produced record
//! *sequence*, so the Rust port is a straightforward lazy iterator (output-identical).

use std::path::Path;
use std::sync::Arc;

use crate::blbutil::{FileIt, Filter, SampleFileIt, VcfFileIt};

use super::{
    BasicGTRec, BitArrayGTRec, GTRec, LowMafDiallelicGTRec, LowMafGTRec, Marker, MarkerParser,
    Samples, VcfHeader, VcfRec, VcfRecGTParser,
};

/// A function mapping `(header, vcf record line, field filter)` to a `GTRec` — the Rust
/// analogue of Java's `TriFunction<VcfHeader, String, MarkerParser, GTRec>`.
pub type RecMapper = fn(&Arc<VcfHeader>, &str, &MarkerParser) -> Box<dyn GTRec>;

/// `VcfIt.TO_LOWMEM_GT_REC` — memory-efficient per-record encoding (low-MAF or bit-packed).
pub fn to_lowmem_gt_rec(h: &Arc<VcfHeader>, s: &str, f: &MarkerParser) -> Box<dyn GTRec> {
    let parser = VcfRecGTParser::new(h.as_ref(), s, f);
    let hlr = parser.hap_list_rep();
    let nonmajor_allele_threshold = hlr.samples().size() >> 7;
    if hlr.nonmajor_allele_cnt() <= nonmajor_allele_threshold {
        if hlr.marker().n_alleles() == 2 {
            Box::new(LowMafDiallelicGTRec::from_hap_list_rep(&hlr))
        } else {
            Box::new(LowMafGTRec::from_hap_list_rep(&hlr))
        }
    } else {
        Box::new(BitArrayGTRec::from_hap_list_rep(&hlr))
    }
}

/// `VcfIt.TO_BASIC_GT_REC` — per-genotype phase status stored.
pub fn to_basic_gt_rec(h: &Arc<VcfHeader>, s: &str, f: &MarkerParser) -> Box<dyn GTRec> {
    Box::new(BasicGTRec::from_parser(&VcfRecGTParser::new(
        h.as_ref(),
        s,
        f,
    )))
}

/// `VcfIt.TO_VCF_REC` — a full `VcfRec`.
pub fn to_vcf_rec(h: &Arc<VcfHeader>, s: &str, f: &MarkerParser) -> Box<dyn GTRec> {
    Box::new(VcfRec::new(h.clone(), s.to_string(), f))
}

/// Port of `vcf/VcfIt.java`.
pub struct VcfIt<I: FileIt<Item = String>> {
    vcf_header: Arc<VcfHeader>,
    it: I,
    mapper: RecMapper,
    field_filter: MarkerParser,
    marker_filter: Filter<Marker>,
    next: Option<String>,
}

/// `VcfIt.head(src, it)` — the leading `#` lines plus the first data line.
pub(crate) fn head<I: Iterator<Item = String>>(src: &str, it: &mut I) -> Vec<String> {
    let mut lines: Vec<String> = Vec::with_capacity(32);
    let mut line = it.next();
    while let Some(l) = line {
        if !l.starts_with('#') {
            line = Some(l);
            break;
        }
        lines.push(l);
        line = it.next();
    }
    match line {
        None => panic!("ERROR: missing VCF data lines ({src})"),
        Some(l) => lines.push(l),
    }
    lines
}

/// `VcfIt.readLine(it)` — the next line, skipping blank lines (returns a trailing blank
/// line at EOF, matching Java).
fn read_line<I: Iterator<Item = String>>(it: &mut I) -> Option<String> {
    let mut line = it.next()?;
    while line.trim().is_empty() {
        match it.next() {
            Some(l) => line = l,
            None => break,
        }
    }
    Some(line)
}

impl<I: FileIt<Item = String>> VcfIt<I> {
    /// `VcfIt.create(strIt, recMapper)` — accept-all sample/marker filters.
    pub fn create(it: I, rec_mapper: RecMapper) -> Self {
        VcfIt::create_filtered(it, Filter::accept_all(), Filter::accept_all(), rec_mapper)
    }

    /// `VcfIt.create(strIt, sampleFilter, markerFilter, recMapper)`.
    pub fn create_filtered(
        mut it: I,
        sample_filter: Filter<String>,
        marker_filter: Filter<Marker>,
        rec_mapper: RecMapper,
    ) -> Self {
        let src = it
            .file()
            .and_then(|p| p.file_name())
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "stdin".to_string());
        let head = head(&src, &mut it);
        let (non_data_lines, first_data) = head.split_at(head.len() - 1);
        let first_data_line = first_data[0].clone();
        let is_diploid = VcfHeader::is_diploid(&first_data_line);
        // storeId=true; do not store qual/filter/info
        let field_filter = MarkerParser::new(true, false, false, false);
        let vcf_header = Arc::new(VcfHeader::new(
            &src,
            non_data_lines,
            &is_diploid,
            &sample_filter,
        ));
        VcfIt {
            vcf_header,
            it,
            mapper: rec_mapper,
            field_filter,
            marker_filter,
            next: Some(first_data_line),
        }
    }
}

impl<I: FileIt<Item = String>> Iterator for VcfIt<I> {
    type Item = Box<dyn GTRec>;

    fn next(&mut self) -> Option<Box<dyn GTRec>> {
        loop {
            let line = self.next.take()?;
            self.next = read_line(&mut self.it);
            let rec = (self.mapper)(&self.vcf_header, &line, &self.field_filter);
            if self.marker_filter.accept(rec.marker()) {
                return Some(rec);
            }
        }
    }
}

impl<I: FileIt<Item = String>> FileIt for VcfIt<I> {
    fn file(&self) -> Option<&Path> {
        self.it.file()
    }
}

impl<I: FileIt<Item = String>> SampleFileIt for VcfIt<I> {
    fn samples(&self) -> &Samples {
        self.vcf_header.samples()
    }
}

impl<I: FileIt<Item = String>> VcfFileIt for VcfIt<I> {
    fn vcf_header(&self) -> &VcfHeader {
        self.vcf_header.as_ref()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blbutil::InputIt;
    use std::io::{BufReader, Cursor};

    const VCF: &str = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\tS1
chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0
chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t1|1\t0|0
chr1\t300\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1
";

    fn line_it(text: &str) -> InputIt {
        InputIt::from_reader(Box::new(BufReader::new(Cursor::new(text.to_owned()))), None)
    }

    #[test]
    fn iterates_records_basic_mapper() {
        let it = VcfIt::create(line_it(VCF), to_basic_gt_rec);
        assert_eq!(it.samples().size(), 2);
        let recs: Vec<Box<dyn GTRec>> = it.collect();
        assert_eq!(recs.len(), 3);
        assert_eq!(recs[0].marker().pos(), 100);
        assert_eq!(recs[1].marker().pos(), 200);
        // record 0: S0=0|1, S1=1|0 -> haps [0,1,1,0]
        assert_eq!(
            (0..4).map(|h| recs[0].get(h)).collect::<Vec<_>>(),
            vec![0, 1, 1, 0]
        );
    }

    #[test]
    fn marker_filter_excludes_records() {
        // exclude marker at chr1:200
        let mf = Filter::predicate(|m: &Marker| m.pos() != 200);
        let it = VcfIt::create_filtered(line_it(VCF), Filter::accept_all(), mf, to_lowmem_gt_rec);
        let positions: Vec<i32> = it.map(|r| r.marker().pos()).collect();
        assert_eq!(positions, vec![100, 300]);
    }

    #[test]
    fn vcf_rec_mapper_and_header() {
        let it = VcfIt::create(line_it(VCF), to_vcf_rec);
        assert_eq!(it.vcf_header().n_samples(), 2);
        let recs: Vec<Box<dyn GTRec>> = it.collect();
        assert_eq!(recs.len(), 3);
        assert_eq!(recs[2].marker().pos(), 300);
    }

    #[test]
    fn skips_blank_lines() {
        let vcf_with_blanks = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0
chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1

chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t1|0
";
        let it = VcfIt::create(line_it(vcf_with_blanks), to_basic_gt_rec);
        let positions: Vec<i32> = it.map(|r| r.marker().pos()).collect();
        assert_eq!(positions, vec![100, 200]);
    }
}
