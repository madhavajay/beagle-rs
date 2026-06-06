//! Port of `vcf/RefIt.java` — an iterator over `RefGTRec`s parsed from a VCF file of
//! phased, non-missing reference genotypes. High-frequency records are re-encoded with the
//! `bref` sequence coder for memory efficiency; the produced record *sequence* (and each
//! record's genotypes) is independent of that encoding.

use std::collections::VecDeque;
use std::path::Path;
use std::rc::Rc;

use crate::blbutil::{BlockLineReader, FileIt, Filter, SampleFileIt, VcfFileIt};
use crate::bref::{SeqCoder3, MAX_NALLELES};

use super::vcf_it::head;
use super::{
    allele_ref_gt_rec_from_parser, Marker, MarkerParser, RefGTRec, Samples, VcfHeader,
    VcfRecGTParser,
};

const DEFAULT_BUFFER_SIZE: i32 = 1 << 10;

/// Port of `vcf/RefIt.java`.
pub struct RefIt<I: FileIt<Item = String>> {
    vcf_header: VcfHeader,
    field_filter: MarkerParser,
    marker_filter: Filter<Marker>,
    reader: BlockLineReader<I>,
    low_freq_buffer: Vec<Option<Box<dyn RefGTRec>>>,
    rec_buffer: VecDeque<Box<dyn RefGTRec>>,
    seq_coder: SeqCoder3,
    max_seq_coded_alleles: i32,
    max_seq_coding_major_cnt: i32,
    last_chrom: i32,
}

fn max_seq_coding_major_cnt(samples: &Samples) -> i32 {
    let n_haps = samples.size() << 1;
    ((n_haps as f32 * crate::bref::COMPRESS_FREQ_THRESHOLD - 1.0_f32) as f64).floor() as i32
}

fn combine(first_data_line: String, lines: Vec<String>) -> Vec<String> {
    let mut v = Vec::with_capacity(lines.len() + 1);
    v.push(first_data_line);
    v.extend(lines);
    v
}

/// `applySeqCoding(RefGTRec rec)`.
fn apply_seq_coding(rec: &dyn RefGTRec, max_seq_coded_alleles: i32, max_major_cnt: i32) -> bool {
    debug_assert!(rec.is_allele_coded());
    if rec.marker().n_alleles() > max_seq_coded_alleles {
        return false;
    }
    let n_haps = rec.size();
    let maj_allele = rec.major_allele();
    let mut maj_cnt = n_haps;
    for a in 0..rec.marker().n_alleles() {
        if a != maj_allele {
            maj_cnt -= rec.allele_count(a);
        }
    }
    maj_cnt <= max_major_cnt
}

impl<I: FileIt<Item = String>> RefIt<I> {
    /// `RefIt.create(FileIt<String> it)`.
    pub fn create(it: I) -> Self {
        RefIt::create_filtered(
            it,
            Filter::accept_all(),
            Filter::accept_all(),
            DEFAULT_BUFFER_SIZE,
        )
    }

    /// `RefIt.create(FileIt<String> it, Filter<String>, Filter<Marker>, int bufferSize)`.
    pub fn create_filtered(
        mut it: I,
        sample_filter: Filter<String>,
        marker_filter: Filter<Marker>,
        block_size: i32,
    ) -> Self {
        assert!(block_size >= 1, "{block_size}");
        let src = it
            .file()
            .and_then(|p| p.file_name())
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "stdin".to_string());
        let head_lines = head(&src, &mut it);
        let (non_data_lines, first_data) = head_lines.split_at(head_lines.len() - 1);
        let first_data_line = first_data[0].clone();
        let is_diploid = VcfHeader::is_diploid(&first_data_line);
        let field_filter = MarkerParser::new(true, false, false, false);
        let vcf_header = VcfHeader::new(&src, non_data_lines, &is_diploid, &sample_filter);
        let seq_coder = SeqCoder3::new(vcf_header.samples().clone());
        let max_seq_coded_alleles = seq_coder.max_n_seq().min(MAX_NALLELES);
        let max_seq_coding_major_cnt = max_seq_coding_major_cnt(vcf_header.samples());
        let reader = BlockLineReader::create(it, block_size, 1);
        let mut ref_it = RefIt {
            vcf_header,
            field_filter,
            marker_filter,
            reader,
            low_freq_buffer: Vec::new(),
            rec_buffer: VecDeque::with_capacity(block_size as usize),
            seq_coder,
            max_seq_coded_alleles,
            max_seq_coding_major_cnt,
            last_chrom: -1,
        };
        ref_it.fill_rec_buffer(Some(first_data_line));
        ref_it
    }

    fn map_line(&self, s: &str) -> Box<dyn RefGTRec> {
        let parser = VcfRecGTParser::new(&self.vcf_header, s, &self.field_filter);
        allele_ref_gt_rec_from_parser(&parser)
    }

    fn parse_lines(&self, lines: &[String]) -> Vec<Box<dyn RefGTRec>> {
        lines
            .iter()
            .map(|s| self.map_line(s))
            .filter(|rec| self.marker_filter.accept(rec.marker()))
            .collect()
    }

    fn fill_rec_buffer(&mut self, mut first_data_line: Option<String>) {
        while self.rec_buffer.is_empty() {
            let mut lines = self.reader.next_block();
            if let Some(fdl) = first_data_line.take() {
                lines = combine(fdl, lines);
            }
            if lines.is_empty() {
                self.flush_compressed_records();
                return;
            }
            let recs = self.parse_lines(&lines);
            for rec in recs {
                let chrom = rec.marker().chrom_index();
                if self.last_chrom == -1 {
                    self.last_chrom = chrom;
                }
                if chrom != self.last_chrom {
                    self.flush_compressed_records();
                    self.last_chrom = chrom;
                }
                if !apply_seq_coding(
                    rec.as_ref(),
                    self.max_seq_coded_alleles,
                    self.max_seq_coding_major_cnt,
                ) {
                    self.low_freq_buffer.push(Some(rec));
                } else {
                    let rec_rc: Rc<dyn RefGTRec> = Rc::from(rec);
                    let mut success = self.seq_coder.add(rec_rc.clone());
                    if !success {
                        self.flush_compressed_records();
                        success = self.seq_coder.add(rec_rc);
                        assert!(success);
                    }
                    self.low_freq_buffer.push(None);
                }
            }
        }
    }

    fn flush_compressed_records(&mut self) {
        let list = self.seq_coder.get_compressed_list();
        let mut list_iter = list.into_iter();
        for slot in self.low_freq_buffer.iter_mut() {
            if slot.is_none() {
                *slot = Some(list_iter.next().expect("compressed record for placeholder"));
            }
        }
        for slot in self.low_freq_buffer.drain(..) {
            self.rec_buffer
                .push_back(slot.expect("filled low-freq slot"));
        }
    }
}

impl<I: FileIt<Item = String>> Iterator for RefIt<I> {
    type Item = Box<dyn RefGTRec>;

    fn next(&mut self) -> Option<Box<dyn RefGTRec>> {
        let first = self.rec_buffer.pop_front()?;
        if self.rec_buffer.is_empty() {
            self.fill_rec_buffer(None);
        }
        Some(first)
    }
}

impl<I: FileIt<Item = String>> FileIt for RefIt<I> {
    fn file(&self) -> Option<&Path> {
        self.reader.file()
    }
}

impl<I: FileIt<Item = String>> SampleFileIt for RefIt<I> {
    fn samples(&self) -> &Samples {
        self.vcf_header.samples()
    }
}

impl<I: FileIt<Item = String>> VcfFileIt for RefIt<I> {
    fn vcf_header(&self) -> &VcfHeader {
        &self.vcf_header
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blbutil::InputIt;
    use std::io::{BufReader, Cursor};

    const VCF: &str = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3
chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0\t0|0\t1|1
chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|1\t0|0
chr1\t300\t.\tA\tT\t.\tPASS\t.\tGT\t1|1\t1|1\t1|1\t0|1
";

    fn ref_it() -> RefIt<InputIt> {
        let it = InputIt::from_reader(Box::new(BufReader::new(Cursor::new(VCF.to_owned()))), None);
        RefIt::create(it)
    }

    // expected per-hap alleles from the genotype strings above
    fn expected(line_alleles: &[(i32, i32)]) -> Vec<i32> {
        line_alleles.iter().flat_map(|&(a, b)| [a, b]).collect()
    }

    #[test]
    fn iterates_reference_records_in_order_with_correct_genotypes() {
        let it = ref_it();
        assert_eq!(it.samples().size(), 4);
        let recs: Vec<Box<dyn RefGTRec>> = it.collect();
        assert_eq!(recs.len(), 3);
        assert_eq!(recs[0].marker().pos(), 100);
        assert_eq!(recs[1].marker().pos(), 200);
        assert_eq!(recs[2].marker().pos(), 300);

        let exp0 = expected(&[(0, 1), (1, 0), (0, 0), (1, 1)]);
        let exp1 = expected(&[(0, 0), (0, 0), (0, 1), (0, 0)]);
        let exp2 = expected(&[(1, 1), (1, 1), (1, 1), (0, 1)]);
        let got = |r: &dyn RefGTRec| (0..8).map(|h| r.get(h)).collect::<Vec<_>>();
        assert_eq!(got(recs[0].as_ref()), exp0);
        assert_eq!(got(recs[1].as_ref()), exp1);
        assert_eq!(got(recs[2].as_ref()), exp2);
        assert!(recs.iter().all(|r| r.is_phased()));
    }

    #[test]
    fn marker_filter_excludes_records() {
        let it = InputIt::from_reader(Box::new(BufReader::new(Cursor::new(VCF.to_owned()))), None);
        let mf = Filter::predicate(|m: &Marker| m.pos() != 200);
        let ri = RefIt::create_filtered(it, Filter::accept_all(), mf, DEFAULT_BUFFER_SIZE);
        let positions: Vec<i32> = ri.map(|r| r.marker().pos()).collect();
        assert_eq!(positions, vec![100, 300]);
    }
}
