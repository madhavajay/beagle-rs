//! Port of `bref/CompressBref3Writer.java` — buffers records and routes each to either
//! sequence-coding (via `SeqCoder3`, for markers with enough non-major-allele carriers) or
//! as-is allele-coding, then emits them in order through an `AsIsBref3Writer`.

use std::io::Write;
use std::rc::Rc;

use crate::vcf::{allele_ref_gt_rec_from_rec, RefGTRec, Samples};

use super::{AsIsBref3Writer, BrefWriter, SeqCoder3, MAX_NALLELES};

/// Port of `bref/CompressBref3Writer.java`.
pub struct CompressBref3Writer<W: Write> {
    max_n_alleles: i32,
    buffer: Vec<Option<Rc<dyn RefGTRec>>>, // `None` marks a sequence-coded record
    seq_coder: SeqCoder3,
    as_is_bref3_writer: AsIsBref3Writer<W>,
    non_maj_threshold: i32,
}

impl<W: Write> CompressBref3Writer<W> {
    /// `new CompressBref3Writer(String program, Samples samples, int maxNSeq, File brefFile)`.
    pub fn new(program: &str, samples: Samples, max_n_seq: i32, writer: W) -> Self {
        CompressBref3Writer {
            max_n_alleles: MAX_NALLELES,
            buffer: Vec::with_capacity(500),
            seq_coder: SeqCoder3::with_max_n_seq(samples.clone(), max_n_seq),
            as_is_bref3_writer: AsIsBref3Writer::new(program, samples, writer),
            non_maj_threshold: (max_n_seq / 4) + 1,
        }
    }

    fn convert_to_seq_coding(&self, rec: &dyn RefGTRec) -> bool {
        // assert rec.is_allele_coded()
        if rec.marker().n_alleles() > self.max_n_alleles {
            return false;
        }
        let maj_allele = rec.major_allele();
        let mut non_maj_cnt = 0;
        for a in 0..rec.marker().n_alleles() {
            if a != maj_allele {
                non_maj_cnt += rec.allele_count(a);
            }
        }
        non_maj_cnt >= self.non_maj_threshold
    }

    fn flush_buffer(&mut self) {
        let list: Vec<Rc<dyn RefGTRec>> = self
            .seq_coder
            .get_compressed_list()
            .into_iter()
            .map(Rc::from)
            .collect();
        let buffer = std::mem::take(&mut self.buffer);
        let mut index = 0;
        for rec_opt in buffer {
            match rec_opt {
                None => {
                    self.as_is_bref3_writer.write(list[index].clone());
                    index += 1;
                }
                Some(rec) => self.as_is_bref3_writer.write(rec),
            }
        }
        debug_assert_eq!(index, list.len());
    }
}

impl<W: Write> BrefWriter for CompressBref3Writer<W> {
    fn samples(&self) -> &Samples {
        self.as_is_bref3_writer.samples()
    }

    fn write(&mut self, rec: Rc<dyn RefGTRec>) {
        let rec: Rc<dyn RefGTRec> = if !rec.is_allele_coded() {
            Rc::from(allele_ref_gt_rec_from_rec(rec.as_ref()))
        } else {
            rec
        };
        if self.convert_to_seq_coding(rec.as_ref()) {
            let success = self.seq_coder.add(rec.clone());
            if !success {
                self.flush_buffer();
                let success = self.seq_coder.add(rec.clone());
                debug_assert!(success);
            }
            self.buffer.push(None);
        } else {
            self.buffer.push(Some(rec));
        }
    }

    fn close(&mut self) {
        self.flush_buffer();
        self.as_is_bref3_writer.close();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bref::Bref3Reader;
    use crate::jdk_io::DataIn;
    use crate::vcf::{
        allele_ref_gt_rec_from_parser, MarkerParser, VcfHeader, VcfRecGTParser, HEADER_PREFIX,
    };
    use std::collections::VecDeque;

    fn header(n: usize) -> VcfHeader {
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..n {
            hdr.push_str(&format!("\tS{s}"));
        }
        VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &vec![true; n],
        )
    }

    fn rec(h: &VcfHeader, line: &str) -> Rc<dyn RefGTRec> {
        let p = VcfRecGTParser::new(h, line, &MarkerParser::new(true, true, true, true));
        Rc::from(allele_ref_gt_rec_from_parser(&p))
    }

    #[test]
    fn write_then_read_round_trips() {
        // 8 samples so that some markers have enough minor-allele carriers to seq-code.
        let h = header(8);
        let lines = [
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\t0|1\t1|0\t0|0\t1|1\t0|1",
            "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0\t0|0\t0|1\t0|0",
            "chr1\t300\t.\tA\tG,T\t.\tPASS\t.\tGT\t0|1\t2|0\t1|1\t0|2\t1|0\t0|0\t2|1\t0|1",
            "chr1\t400\t.\tC\tA\t.\tPASS\t.\tGT\t1|1\t1|0\t0|1\t1|1\t0|0\t1|0\t0|1\t1|1",
        ];
        let recs: Vec<Rc<dyn RefGTRec>> = lines.iter().map(|l| rec(&h, l)).collect();
        let samples = recs[0].samples().clone();

        let mut buf: Vec<u8> = Vec::new();
        {
            let mut w = CompressBref3Writer::new("prog", samples.clone(), 7, &mut buf);
            for r in &recs {
                w.write(r.clone());
            }
            w.close();
        }
        assert!(!buf.is_empty());

        // read it back and confirm the genotypes match
        let mut di = DataIn::new(std::io::Cursor::new(buf));
        let mut reader = Bref3Reader::new(None, &mut di);
        assert_eq!(reader.samples(), &samples);
        let mut buffer: VecDeque<Box<dyn RefGTRec>> = VecDeque::new();
        reader.read_block(&mut di, &mut buffer);
        let mut read: Vec<Box<dyn RefGTRec>> = Vec::new();
        while let Some(rec) = buffer.pop_front() {
            read.push(rec);
            if buffer.is_empty() {
                reader.read_block(&mut di, &mut buffer);
            }
        }
        assert_eq!(read.len(), recs.len());
        for (orig, got) in recs.iter().zip(read.iter()) {
            assert_eq!(got.marker(), orig.marker());
            for hap in 0..samples.size() * 2 {
                assert_eq!(
                    got.get(hap),
                    orig.get(hap),
                    "marker {}",
                    orig.marker().pos()
                );
            }
        }
    }
}
