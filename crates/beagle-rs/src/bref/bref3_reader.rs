//! Port of `bref/Bref3Reader.java` — reads bref3 data blocks into `RefGTRec`s.
//!
//! A block is: `int nRecs`, `UTF chrom`, `unsignedShort nSeq`, a hap→seq byte map for all
//! unfiltered haplotypes, then `nRecs` records. Each record is a marker, a `byte` flag
//! (0 = sequence-coded `HapRefGTRec`, 1 = allele-coded), and the coded genotype data.

use std::collections::VecDeque;
use std::io::{self, Read};
use std::path::Path;
use std::sync::OnceLock;

use crate::beagleutil::ChromIds;
use crate::blbutil::{consts, Filter, Utilities};
use crate::ints::{CharArray, UnsignedByteArray};
use crate::jdk_io::DataIn;
use crate::vcf::{
    allele_ref_gt_rec_from_components, HapRefGTRec, Marker, MarkerParser, RefGTRec, Samples,
};

use super::{read_string_array, Bref3Header};

const READ_ERR: &str = "Error reading file";

/// Port of `bref/Bref3Reader.java`.
pub struct Bref3Reader {
    marker_filter: Filter<Marker>,
    program: String,
    included_hap_indices: Vec<i32>,
    inv_included_hap_indices: Vec<i32>,
    samples: Samples,
    byte_buffer: Vec<u8>,
}

impl Bref3Reader {
    /// `new Bref3Reader(File source, DataInput dataIn)` — accept-all filters.
    pub fn new<R: Read>(source: Option<&Path>, di: &mut DataIn<R>) -> Self {
        Bref3Reader::with_filters(source, di, Filter::accept_all(), Filter::accept_all())
    }

    /// `new Bref3Reader(File, DataInput, Filter<String>, Filter<Marker>)`.
    pub fn with_filters<R: Read>(
        source: Option<&Path>,
        di: &mut DataIn<R>,
        sample_filter: Filter<String>,
        marker_filter: Filter<Marker>,
    ) -> Self {
        let bref_header = Bref3Header::new(source, di, &sample_filter);
        let program = bref_header.program().to_string();
        let included_hap_indices = bref_header.filtered_hap_indices();
        let inv_included_hap_indices = bref_header.inv_filtered_hap_indices();
        let samples = bref_header.samples().clone();
        let byte_buffer = vec![0u8; 2 * inv_included_hap_indices.len()];
        Bref3Reader {
            marker_filter,
            program,
            included_hap_indices,
            inv_included_hap_indices,
            samples,
            byte_buffer,
        }
    }

    /// `samples()`.
    pub fn samples(&self) -> &Samples {
        &self.samples
    }

    /// `program()`.
    pub fn program(&self) -> &str {
        &self.program
    }

    /// `readBlock(DataInput bref, Collection<RefGTRec> buffer)` — reads blocks until the
    /// buffer is non-empty or the end-of-data sentinel (`nRecs == 0`) is reached.
    pub fn read_block<R: Read>(
        &mut self,
        di: &mut DataIn<R>,
        buffer: &mut VecDeque<Box<dyn RefGTRec>>,
    ) {
        if let Err(e) = self.read_block_impl(di, buffer) {
            Utilities::exit(&format!("{READ_ERR}: {e}"));
        }
    }

    fn read_block_impl<R: Read>(
        &mut self,
        di: &mut DataIn<R>,
        buffer: &mut VecDeque<Box<dyn RefGTRec>>,
    ) -> io::Result<()> {
        let mut n_recs = i32::MAX;
        while buffer.is_empty() && n_recs != 0 {
            n_recs = di.read_int()?;
            if n_recs != 0 {
                self.read_block_n(di, buffer, n_recs)?;
            }
        }
        Ok(())
    }

    fn read_block_n<R: Read>(
        &mut self,
        di: &mut DataIn<R>,
        buffer: &mut VecDeque<Box<dyn RefGTRec>>,
        n_recs: i32,
    ) -> io::Result<()> {
        let chrom = di.read_utf()?;
        let chrom_index = ChromIds::instance().get_index(&chrom);
        let n_seq = di.read_unsigned_short()?;
        di.read_fully(&mut self.byte_buffer)?;
        let hap_to_seq = self.hap_to_seq();
        for _ in 0..n_recs {
            let marker = read_marker(di, chrom_index)?;
            let flag = di.read_byte()?;
            let rec: Box<dyn RefGTRec> = match flag {
                0 => self.read_hap_record(di, marker, &hap_to_seq, n_seq)?,
                1 => self.read_allele_record(di, marker)?,
                _ => Utilities::exit(READ_ERR),
            };
            if self.marker_filter.accept(rec.marker()) {
                buffer.push_back(rec);
            }
        }
        Ok(())
    }

    fn hap_to_seq(&self) -> CharArray {
        let values: Vec<i32> = self
            .included_hap_indices
            .iter()
            .map(|&hap| {
                let offset = (hap << 1) as usize;
                let b1 = self.byte_buffer[offset] as i32;
                let b2 = self.byte_buffer[offset + 1] as i32;
                (b1 << 8) + b2
            })
            .collect();
        CharArray::from_ints(&values)
    }

    fn read_hap_record<R: Read>(
        &mut self,
        di: &mut DataIn<R>,
        marker: Marker,
        hap_to_seq: &CharArray,
        n_seq: i32,
    ) -> io::Result<Box<dyn RefGTRec>> {
        di.read_fully(&mut self.byte_buffer[0..n_seq as usize])?;
        let seq_to_allele = UnsignedByteArray::from_bytes_range(&self.byte_buffer, 0, n_seq);
        Ok(Box::new(HapRefGTRec::new(
            marker,
            self.samples.clone(),
            hap_to_seq,
            &seq_to_allele,
        )))
    }

    fn read_allele_record<R: Read>(
        &self,
        di: &mut DataIn<R>,
        marker: Marker,
    ) -> io::Result<Box<dyn RefGTRec>> {
        let n_alleles = marker.n_alleles();
        let mut hap_indices: Vec<Option<Vec<i32>>> = Vec::with_capacity(n_alleles as usize);
        for _ in 0..n_alleles {
            hap_indices.push(self.read_allele_coded_hap_list(di)?);
        }
        Ok(allele_ref_gt_rec_from_components(
            marker,
            self.samples.clone(),
            hap_indices,
        ))
    }

    fn read_allele_coded_hap_list<R: Read>(
        &self,
        di: &mut DataIn<R>,
    ) -> io::Result<Option<Vec<i32>>> {
        let length = di.read_int()?;
        if length == -1 {
            Ok(None)
        } else {
            let mut ia = Vec::with_capacity(length as usize);
            for _ in 0..length {
                let hap = di.read_int()?;
                if self.inv_included_hap_indices[hap as usize] >= 0 {
                    ia.push(self.inv_included_hap_indices[hap as usize]);
                }
            }
            Ok(Some(ia))
        }
    }
}

/// `Bref3Reader.readMarker(DataInput, int chromIndex)`.
fn read_marker<R: Read>(di: &mut DataIn<R>, chrom_index: i32) -> io::Result<Marker> {
    let pos = di.read_int()?;
    let id = read_byte_length_string_array_and_join(di, consts::SEMICOLON)?;
    let allele_code = di.read_byte()? as i32;
    let (ref_and_alt, end) = if allele_code == -1 {
        let str_alleles = read_string_array(di)?.expect("allele array");
        let end = di.read_int()?;
        (ref_and_alt_fields(&str_alleles), end)
    } else {
        let n_alleles = 1 + (allele_code & 0b11);
        let perm_index = allele_code >> 2;
        let str_alleles = allele_string(perm_index, n_alleles);
        (ref_and_alt_fields(&str_alleles), -1)
    };
    let vcf_rec_prefix = vcf_rec_prefix(chrom_index, pos, &id, &ref_and_alt, end);
    Ok(Marker::instance(
        &vcf_rec_prefix,
        &MarkerParser::new(true, false, false, false),
    ))
}

fn ref_and_alt_fields(str_alleles: &[String]) -> String {
    let mut sb = String::new();
    sb.push_str(&str_alleles[0]);
    if str_alleles.len() == 1 {
        sb.push(consts::TAB);
        sb.push(consts::MISSING_DATA_CHAR);
    } else {
        for (j, allele) in str_alleles.iter().enumerate().skip(1) {
            sb.push(if j == 1 { consts::TAB } else { consts::COMMA });
            sb.push_str(allele);
        }
    }
    sb
}

fn vcf_rec_prefix(chrom: i32, pos: i32, id: &str, ref_and_alt_fields: &str, end: i32) -> String {
    let mut sb = String::with_capacity(64);
    sb.push_str(&ChromIds::instance().id(chrom));
    sb.push(consts::TAB);
    sb.push_str(&pos.to_string());
    sb.push(consts::TAB);
    sb.push_str(id);
    sb.push(consts::TAB);
    sb.push_str(ref_and_alt_fields);
    sb.push(consts::TAB);
    sb.push(consts::MISSING_DATA_CHAR); // QUAL
    sb.push(consts::TAB);
    sb.push(consts::MISSING_DATA_CHAR); // FILTER
    sb.push(consts::TAB);
    if end >= 0 {
        sb.push_str("END="); // INFO
        sb.push_str(&end.to_string());
    } else {
        sb.push(consts::MISSING_DATA_CHAR); // INFO
    }
    sb.push(consts::TAB);
    sb.push_str("GT");
    sb
}

fn read_byte_length_string_array_and_join<R: Read>(
    di: &mut DataIn<R>,
    delim: char,
) -> io::Result<String> {
    let length = di.read_unsigned_byte()?;
    if length <= 0 {
        Ok(consts::MISSING_DATA_STRING.to_string())
    } else {
        let mut sb = String::new();
        for j in 0..length {
            if j > 0 {
                sb.push(delim);
            }
            sb.push_str(&di.read_utf()?);
        }
        Ok(sb)
    }
}

/// The 24 permutations of `["A","C","G","T"]` in lexicographic order, generated by the same
/// recursion as Java's `Bref3Reader.permute`.
pub(crate) fn snv_perms() -> &'static Vec<Vec<&'static str>> {
    static PERMS: OnceLock<Vec<Vec<&'static str>>> = OnceLock::new();
    PERMS.get_or_init(|| {
        let mut perms = Vec::with_capacity(24);
        permute(Vec::new(), &["A", "C", "G", "T"], &mut perms);
        perms
    })
}

fn permute(start: Vec<&'static str>, end: &[&'static str], perms: &mut Vec<Vec<&'static str>>) {
    if end.is_empty() {
        perms.push(start);
    } else {
        for j in 0..end.len() {
            let mut new_start = start.clone();
            new_start.push(end[j]);
            let mut new_end = Vec::with_capacity(end.len() - 1);
            new_end.extend_from_slice(&end[..j]);
            new_end.extend_from_slice(&end[j + 1..]);
            permute(new_start, &new_end, perms);
        }
    }
}

/// `alleleString(int permIndex, int length)` — first `length` bases of permutation
/// `permIndex`.
fn allele_string(perm_index: i32, length: i32) -> Vec<String> {
    snv_perms()[perm_index as usize][..length as usize]
        .iter()
        .map(|s| s.to_string())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jdk_io::DataOut;
    use std::io::Cursor;

    #[test]
    fn snv_perms_lexicographic_order() {
        let perms = snv_perms();
        assert_eq!(perms.len(), 24);
        assert_eq!(perms[0], vec!["A", "C", "G", "T"]);
        assert_eq!(perms[1], vec!["A", "C", "T", "G"]);
        assert_eq!(perms[23], vec!["T", "G", "C", "A"]);
    }

    #[test]
    fn ref_and_alt_formatting() {
        assert_eq!(
            ref_and_alt_fields(&["A".to_string(), "C".to_string()]),
            "A\tC"
        );
        assert_eq!(
            ref_and_alt_fields(&["A".to_string(), "C".to_string(), "G".to_string()]),
            "A\tC,G"
        );
        // single allele -> ALT is missing
        assert_eq!(ref_and_alt_fields(&["A".to_string()]), "A\t.");
    }

    #[test]
    fn reads_allele_coded_snv_marker() {
        let chrom_index = ChromIds::instance().get_index("chrBREFR");
        // pos, id length 0 (missing), allele_code = (permIndex 0 << 2) | (nAlleles 2 - 1) = 1
        let mut out = DataOut::new(Vec::new());
        out.write_int(5000).unwrap();
        out.write_byte(0).unwrap(); // id: length 0 -> "."
        out.write_byte(1).unwrap(); // allele_code 1 -> biallelic, perm 0 = A,C
        let mut di = DataIn::new(Cursor::new(out.into_inner()));
        let m = read_marker(&mut di, chrom_index).unwrap();
        assert_eq!(m.chrom(), "chrBREFR");
        assert_eq!(m.pos(), 5000);
        assert_eq!(m.n_alleles(), 2);
        assert_eq!(crate::vcf::marker_utils::alleles(&m), vec!["A", "C"]);
    }

    #[test]
    fn reads_non_snv_marker_with_end() {
        let chrom_index = ChromIds::instance().get_index("chrBREFR2");
        let mut out = DataOut::new(Vec::new());
        out.write_int(7000).unwrap();
        // id: 1 string "rsX"
        out.write_byte(1).unwrap();
        out.write_utf("rsX").unwrap();
        out.write_byte(0xFF).unwrap(); // allele_code -1 -> general alleles follow
        out.write_int(2).unwrap(); // 2 alleles
        out.write_utf("AT").unwrap();
        out.write_utf("A").unwrap();
        out.write_int(7100).unwrap(); // END
        let mut di = DataIn::new(Cursor::new(out.into_inner()));
        let m = read_marker(&mut di, chrom_index).unwrap();
        assert_eq!(m.pos(), 7000);
        assert_eq!(m.id(), "rsX");
        assert_eq!(crate::vcf::marker_utils::alleles(&m), vec!["AT", "A"]);
        assert_eq!(m.end_value(), "7100");
    }
}
