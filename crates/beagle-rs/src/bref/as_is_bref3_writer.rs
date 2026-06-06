//! Port of `bref/AsIsBref3Writer.java` — writes phased, non-missing genotypes to a bref3 file,
//! preserving each record's allele-coded vs sequence-coded representation. Byte-exact: the
//! running `bytes_written` count drives the block index offsets written in the file footer.

use std::io::Write;
use std::rc::Rc;

use crate::beagleutil::ChromIds;
use crate::blbutil::{consts, Utilities};
use crate::ints::IntList;
use crate::jdk_io::{utf_byte_len, DataOut};
use crate::vcf::{marker_utils, Marker, RefGTRec, Samples};

use super::bref3_reader::snv_perms;
use super::{
    BrefBlock, BrefWriter, ALLELE_CODED, END_OF_DATA, END_OF_INDEX, MAGIC_NUMBER_V3, SEQ_CODED,
};

const MAX_SAMPLES: i32 = (1 << 29) - 1;
const WRITE_ERR: &str = "Error writing file";
const CONTIGUITY_ERR: &str = "Error: chromosomes not contiguous";

/// Port of `bref/AsIsBref3Writer.java`.
pub struct AsIsBref3Writer<W: Write> {
    bref_out: DataOut<W>,
    bytes_written: i64,
    samples: Samples,
    n_haps: i32,
    rec_buffer: Vec<Rc<dyn RefGTRec>>,
    index: Vec<BrefBlock>,
    last_chrom_index: i32,
    hap2_seq: Option<usize>,
}

fn is_base(s: &str) -> bool {
    matches!(s, "A" | "C" | "G" | "T")
}

fn is_snv(allele_list: &[String]) -> bool {
    allele_list.iter().all(|a| is_base(a))
}

/// `ALLELES_COMP` — compares two allele lists by the first character of each element, then by
/// length (shorter first) when one is a prefix of the other.
fn alleles_cmp(o1: &[&str], o2: &[String]) -> i32 {
    let n = o1.len().min(o2.len());
    for k in 0..n {
        let c1 = o1[k].as_bytes()[0];
        let c2 = o2[k].as_bytes()[0];
        if c1 != c2 {
            return if c1 < c2 { -1 } else { 1 };
        }
    }
    match o1.len().cmp(&o2.len()) {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Greater => 1,
        std::cmp::Ordering::Equal => 0,
    }
}

/// `Arrays.binarySearch(SNV_PERMS, alleles, ALLELES_COMP)` semantics.
fn binary_search_perms(perms: &[Vec<&str>], key: &[String]) -> i32 {
    let mut lo = 0i32;
    let mut hi = perms.len() as i32 - 1;
    while lo <= hi {
        let mid = (lo + hi) >> 1;
        let cmp = alleles_cmp(&perms[mid as usize], key);
        match cmp.cmp(&0) {
            std::cmp::Ordering::Less => lo = mid + 1,
            std::cmp::Ordering::Greater => hi = mid - 1,
            std::cmp::Ordering::Equal => return mid,
        }
    }
    -(lo + 1)
}

fn snv_code(alleles: &[String]) -> i32 {
    let x = binary_search_perms(snv_perms(), alleles);
    let x = if x < 0 { -x - 1 } else { x };
    (x << 2) + (alleles.len() as i32 - 1)
}

fn extract_end(marker: &Marker) -> i32 {
    let info = marker.info();
    let mut index: usize = 4; // start of base coordinate if info.startsWith("END=")
    if !info.starts_with("END=") {
        match info.find(";END=") {
            Some(p) => index = p + 5,
            None => return -1,
        }
    }
    let end_index = info[index..]
        .find(consts::SEMICOLON)
        .map(|p| index + p)
        .unwrap_or(info.len());
    info[index..end_index].parse::<i32>().expect("END value")
}

impl<W: Write> AsIsBref3Writer<W> {
    /// `new AsIsBref3Writer(String program, Samples samples, File brefFile)` — `writer` is the
    /// output sink (a file/stdout in production, any `Write` for testing).
    pub fn new(program: &str, samples: Samples, writer: W) -> Self {
        if samples.size() > MAX_SAMPLES {
            panic!("{}", samples.size());
        }
        let n_haps = 2 * samples.size();
        let mut w = AsIsBref3Writer {
            bref_out: DataOut::new(writer),
            bytes_written: 0,
            samples,
            n_haps,
            rec_buffer: Vec::with_capacity(500),
            index: Vec::with_capacity(500),
            last_chrom_index: -1,
            hap2_seq: None,
        };
        if w.bref_out.write_int(MAGIC_NUMBER_V3).is_err() {
            Utilities::exit(WRITE_ERR);
        }
        w.bytes_written += 4;
        let ids = w.samples.ids();
        w.write_string(program);
        w.write_string_array(&ids);
        w
    }

    fn write_string(&mut self, s: &str) {
        if self.bref_out.write_utf(s).is_err() {
            Utilities::exit(WRITE_ERR);
        }
        self.bytes_written += utf_byte_len(s) as i64;
    }

    fn write_string_array(&mut self, sa: &[String]) {
        if self.bref_out.write_int(sa.len() as i32).is_err() {
            Utilities::exit(WRITE_ERR);
        }
        self.bytes_written += 4;
        for s in sa {
            self.write_string(s);
        }
    }

    fn start_new_block(&mut self, rec: &Rc<dyn RefGTRec>) -> bool {
        let mut start_new_block = false;
        let chrom_index = rec.marker().chrom_index();
        if chrom_index != self.last_chrom_index {
            self.last_chrom_index = chrom_index;
            self.hap2_seq = None;
            start_new_block = true;
        }
        if !rec.is_allele_coded() {
            let key = rec.seq_block_key();
            if self.hap2_seq.is_none() {
                self.hap2_seq = key;
            } else if key != self.hap2_seq {
                self.hap2_seq = key;
                start_new_block = true;
            }
        }
        start_new_block
    }

    fn write_and_clear_buffer(&mut self) {
        if self.rec_buffer.is_empty() {
            return;
        }
        let buffer = std::mem::take(&mut self.rec_buffer);
        let m = buffer[0].marker().clone();
        self.index
            .push(BrefBlock::new(m.chrom_index(), m.pos(), self.bytes_written));
        if self.bref_out.write_int(buffer.len() as i32).is_err() {
            Utilities::exit(WRITE_ERR);
        }
        self.bytes_written += 4;
        self.write_string(&m.chrom());
        self.write_hap_to_seq(&buffer);
        for rec in &buffer {
            if rec.is_allele_coded() {
                self.write_allele_coded_rec(rec.as_ref());
            } else {
                self.write_seq_coded_rec(rec.as_ref());
            }
        }
    }

    fn write_hap_to_seq(&mut self, buffer: &[Rc<dyn RefGTRec>]) {
        let rec = buffer.iter().find(|r| !r.is_allele_coded());
        match rec {
            None => {
                let _ = self.bref_out.write_char(0);
                for _ in 0..self.n_haps {
                    let _ = self.bref_out.write_char(0);
                }
            }
            Some(rec) => {
                let hap_to_seq = rec.map(0);
                let seq_to_allele = rec.map(1);
                let _ = self.bref_out.write_char(seq_to_allele.size());
                for j in 0..hap_to_seq.size() {
                    let _ = self.bref_out.write_char(hap_to_seq.get(j));
                }
            }
        }
        self.bytes_written += (self.n_haps as i64 + 1) * 2; // Character.BYTES
    }

    fn write_seq_coded_rec(&mut self, rec: &dyn RefGTRec) {
        if rec.marker().n_alleles() >= 256 {
            Utilities::exit("ERROR: more than 256 alleles");
        }
        let seq2_allele = rec.map(1);
        self.write_marker(rec.marker());
        let _ = self.bref_out.write_byte(SEQ_CODED as i32);
        for j in 0..seq2_allele.size() {
            let _ = self.bref_out.write_byte(seq2_allele.get(j));
        }
        self.bytes_written += seq2_allele.size() as i64 + 1; // (size + 1) * Byte.BYTES(1)
    }

    fn write_allele_coded_rec(&mut self, rec: &dyn RefGTRec) {
        let n_alleles = rec.marker().n_alleles();
        let major_allele = rec.major_allele();
        self.write_marker(rec.marker());
        let _ = self.bref_out.write_byte(ALLELE_CODED as i32);
        self.bytes_written += 1;
        for a in 0..n_alleles {
            if a == major_allele {
                let _ = self.bref_out.write_int(-1);
                self.bytes_written += 4;
            } else {
                let al_cnt = rec.allele_count(a);
                let _ = self.bref_out.write_int(al_cnt);
                for c in 0..al_cnt {
                    let _ = self.bref_out.write_int(rec.hap_index(a, c));
                }
                self.bytes_written += (al_cnt as i64 + 1) * 4;
            }
        }
    }

    fn write_marker(&mut self, marker: &Marker) {
        let ids = marker_utils::ids(marker);
        let n_ids = ids.len().min(255);
        let _ = self.bref_out.write_int(marker.pos());
        let _ = self.bref_out.write_byte(n_ids as i32);
        self.bytes_written += 4 + 1;
        for id in ids.iter().take(n_ids) {
            self.write_string(id);
        }
        let allele_list = marker_utils::alleles(marker);
        let allele_code = if is_snv(&allele_list) {
            snv_code(&allele_list)
        } else {
            -1
        };
        let _ = self.bref_out.write_byte(allele_code);
        self.bytes_written += 1;
        if allele_code == -1 {
            self.write_string_array(&allele_list);
            let _ = self.bref_out.write_int(extract_end(marker));
            self.bytes_written += 4;
        }
    }

    fn write_index(&mut self) {
        let index = std::mem::take(&mut self.index);
        self.write_index_chroms(&index);
        let mut last_chr_index = -1;
        for bb in &index {
            let mut offset = bb.offset();
            let ci = bb.chrom_index();
            if ci != last_chr_index {
                last_chr_index = ci;
                offset = -offset;
            }
            let _ = self.bref_out.write_long(offset);
            let _ = self.bref_out.write_int(bb.pos());
        }
        let _ = self.bref_out.write_long(END_OF_INDEX);
        self.bytes_written += (8 + 4) * index.len() as i64 + 8;
    }

    fn write_index_chroms(&mut self, index: &[BrefBlock]) {
        let mut last_chr_index = -1;
        let mut chrom_list: Vec<String> = Vec::new();
        let mut first_chrom_block = IntList::new();
        for (j, bb) in index.iter().enumerate() {
            let chrom_index = bb.chrom_index();
            if chrom_index != last_chr_index {
                chrom_list.push(ChromIds::instance().id(chrom_index));
                first_chrom_block.add(j as i32);
                last_chr_index = chrom_index;
            }
        }
        let unique: std::collections::HashSet<&String> = chrom_list.iter().collect();
        if chrom_list.len() != unique.len() {
            Utilities::exit(CONTIGUITY_ERR);
        }
        self.write_string_array(&chrom_list);
        for j in 0..first_chrom_block.size() {
            let _ = self.bref_out.write_int(first_chrom_block.get(j));
        }
    }
}

impl<W: Write> BrefWriter for AsIsBref3Writer<W> {
    fn samples(&self) -> &Samples {
        &self.samples
    }

    fn write(&mut self, rec: Rc<dyn RefGTRec>) {
        if rec.samples() != &self.samples {
            Utilities::exit("ERROR: inconsistent data");
        }
        if self.start_new_block(&rec) {
            self.write_and_clear_buffer();
        }
        self.rec_buffer.push(rec);
    }

    fn close(&mut self) {
        self.write_and_clear_buffer();
        let _ = self.bref_out.write_int(END_OF_DATA);
        self.bytes_written += 4;
        let index_offset = self.bytes_written;
        self.write_index();
        let _ = self.bref_out.write_long(index_offset);
        self.bytes_written += 8;
        if self.bref_out.flush().is_err() {
            Utilities::exit("Error closing file");
        }
    }
}
