//! Port of `vcf/Marker.java` — a VCF record's CHROM/POS/ID/REF/ALT/QUAL/FILTER/INFO,
//! stored compactly: `chrom_index`, `pos`, a packed `field_info` bit flag, and a packed
//! `fields` string. `field_info` is `u16` here (Java `short`) so bit tests are clean.
//!
//! Deferred: `writeNonPosFields`/`readNonPosFields` (bref3 serialization — ported with
//! the `bref` package) and the private "UNDER CONSTRUCTION" `targToRefAllele` family
//! (dead code with no observable behavior).

use crate::beagleutil::ChromIds;
use crate::blbutil::consts;
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};

use super::marker_utils;
use super::MarkerParser;

pub(crate) const STORED_N_ALLELES_MASK: u16 = 0xff;
const INDEXED_N_ALLELES_MASK: u16 = 0b111;
const SNV_INDEX_MASK: u16 = 0x7f;

pub(crate) const ID_STORED: u16 = 1 << 15;
pub(crate) const ALLELES_STORED: u16 = 1 << 14;
pub(crate) const QUAL_STORED: u16 = 1 << 13;
pub(crate) const FILTER_STORED: u16 = 1 << 12;
pub(crate) const INFO_STORED: u16 = 1 << 11;
pub(crate) const END_STORED: u16 = 1 << 10;
/// All "stored" flag bits — set iff the `fields` string is non-empty (used by the
/// deferred bref3 serialization).
#[allow(dead_code)]
pub(crate) const FLAGS: u16 =
    ID_STORED | ALLELES_STORED | QUAL_STORED | FILTER_STORED | INFO_STORED | END_STORED;

/// `String.indexOf(ch, from)` — byte index or -1.
fn index_of(s: &str, ch: char, from: usize) -> i32 {
    s[from..].find(ch).map_or(-1, |i| (from + i) as i32)
}

/// Port of `vcf/Marker.java`.
#[derive(Clone, Debug)]
pub struct Marker {
    chrom_index: i16,
    pos: i32,
    field_info: u16,
    fields: Option<String>,
}

impl Marker {
    pub(crate) fn from_parts(
        chrom_index: i16,
        pos: i32,
        field_info: u16,
        fields: Option<String>,
    ) -> Self {
        Marker {
            chrom_index,
            pos,
            field_info,
            fields,
        }
    }

    /// `Marker.instance(String rec, MarkerParser markerParser)`.
    pub fn instance(rec: &str, marker_parser: &MarkerParser) -> Marker {
        let tabs = marker_utils::first_8_tab_indices(rec);
        let chrom_index = marker_utils::chrom_index(rec, &rec[0..tabs[0] as usize]);
        let pos: i32 = rec[(tabs[0] + 1) as usize..tabs[1] as usize]
            .parse()
            .expect("VCF POS is a parsable integer");
        let mut sb = String::new();
        let field_info = marker_parser.store_marker_fields(rec, &mut sb, &tabs);
        let fields = if sb.is_empty() { None } else { Some(sb) };
        Marker::from_parts(chrom_index, pos, field_info, fields)
    }

    fn fields(&self) -> &str {
        self.fields.as_deref().unwrap_or("")
    }

    /// `chrom()`.
    pub fn chrom(&self) -> String {
        ChromIds::instance().id(self.chrom_index as i32)
    }

    /// `chromIndex()`.
    pub fn chrom_index(&self) -> i32 {
        self.chrom_index as i32
    }

    /// `pos()`.
    pub fn pos(&self) -> i32 {
        self.pos
    }

    /// `hasIdData()`.
    pub fn has_id_data(&self) -> bool {
        self.field_info & ID_STORED == ID_STORED
    }

    /// `hasQualData()`.
    pub fn has_qual_data(&self) -> bool {
        self.field_info & QUAL_STORED == QUAL_STORED
    }

    /// `hasFilterData()`.
    pub fn has_filter_data(&self) -> bool {
        self.field_info & FILTER_STORED == FILTER_STORED
    }

    /// `hasInfoData()`.
    pub fn has_info_data(&self) -> bool {
        self.field_info & INFO_STORED == INFO_STORED
    }

    /// `hasEndValue()`.
    pub fn has_end_value(&self) -> bool {
        self.field_info & END_STORED == END_STORED
    }

    /// `id()`.
    pub fn id(&self) -> String {
        if self.has_id_data() {
            let fields = self.fields();
            let end = index_of(fields, consts::TAB, 0);
            if end < 0 {
                fields.to_string()
            } else {
                fields[0..end as usize].to_string()
            }
        } else {
            consts::MISSING_DATA_STRING.to_string()
        }
    }

    /// `alleles()` — tab-separated REF and ALT fields.
    pub fn alleles(&self) -> String {
        if self.field_info & ALLELES_STORED == ALLELES_STORED {
            let fields = self.fields();
            let mut start = 0usize;
            if self.has_id_data() {
                start = (index_of(fields, consts::TAB, 0) + 1) as usize;
            }
            let end_ref = index_of(fields, consts::TAB, start);
            let end_alt = index_of(fields, consts::TAB, (end_ref + 1) as usize);
            if end_alt < 0 {
                fields[start..].to_string()
            } else {
                fields[start..end_alt as usize].to_string()
            }
        } else {
            let snv_index = ((self.field_info >> 3) & SNV_INDEX_MASK) as usize;
            let n_alleles = (self.field_info & INDEXED_N_ALLELES_MASK) as i32;
            let perms = marker_utils::snv_perms();
            if n_alleles == 1 {
                perms[snv_index].clone()
            } else {
                let length = ((n_alleles << 1) - 1) as usize;
                perms[snv_index][0..length].to_string()
            }
        }
    }

    /// `nRefBases()`.
    pub fn n_ref_bases(&self) -> i32 {
        if self.field_info & ALLELES_STORED == ALLELES_STORED {
            let fields = self.fields();
            let mut start = 0i32;
            if self.has_id_data() {
                start = index_of(fields, consts::TAB, 0) + 1;
            }
            let end = index_of(fields, consts::TAB, start as usize);
            end - start
        } else {
            1
        }
    }

    /// `nAlleles()`.
    pub fn n_alleles(&self) -> i32 {
        if self.field_info & ALLELES_STORED == ALLELES_STORED {
            (self.field_info & STORED_N_ALLELES_MASK) as i32
        } else {
            (self.field_info & INDEXED_N_ALLELES_MASK) as i32
        }
    }

    /// `bitsPerAllele()`.
    pub fn bits_per_allele(&self) -> i32 {
        (32 - ((self.n_alleles() - 1) as u32).leading_zeros()) as i32
    }

    fn qual_start_index(&self) -> i32 {
        let fields = self.fields();
        let mut start = 0i32;
        if self.has_id_data() {
            start = index_of(fields, consts::TAB, 0) + 1;
        }
        if self.field_info & ALLELES_STORED == ALLELES_STORED {
            start = index_of(fields, consts::TAB, start as usize) + 1;
            start = index_of(fields, consts::TAB, start as usize) + 1;
        }
        start
    }

    fn extract_field(&self, start: i32, end_delim: char) -> String {
        let fields = self.fields();
        let end = index_of(fields, end_delim, start as usize);
        if end < 0 {
            fields[start as usize..].to_string()
        } else {
            fields[start as usize..end as usize].to_string()
        }
    }

    /// `qual()`.
    pub fn qual(&self) -> String {
        if self.has_qual_data() {
            self.extract_field(self.qual_start_index(), consts::TAB)
        } else {
            consts::MISSING_DATA_STRING.to_string()
        }
    }

    /// `filter()`.
    pub fn filter(&self) -> String {
        if self.has_filter_data() {
            let mut start = self.qual_start_index();
            if self.has_qual_data() {
                start = index_of(self.fields(), consts::TAB, start as usize) + 1;
            }
            self.extract_field(start, consts::TAB)
        } else {
            consts::MISSING_DATA_STRING.to_string()
        }
    }

    /// `info()`.
    pub fn info(&self) -> String {
        if self.has_info_data() || self.has_end_value() {
            let mut start = self.qual_start_index();
            if self.has_qual_data() {
                start = index_of(self.fields(), consts::TAB, start as usize) + 1;
            }
            if self.has_filter_data() {
                start = index_of(self.fields(), consts::TAB, start as usize) + 1;
            }
            self.extract_field(start, consts::TAB)
        } else {
            consts::MISSING_DATA_STRING.to_string()
        }
    }

    /// `endValue()` — the INFO/END value, or `""` if undefined.
    pub fn end_value(&self) -> String {
        if self.has_end_value() {
            let fields = self.fields();
            let mut start = self.qual_start_index();
            if self.has_qual_data() {
                start = index_of(fields, consts::TAB, start as usize) + 1;
            }
            if self.has_filter_data() {
                start = index_of(fields, consts::TAB, start as usize) + 1;
            }
            // Java: fields.indexOf("END=", startValue); assert startKey >= 0
            let start_key = fields[start as usize..]
                .find("END=")
                .map(|i| start as usize + i)
                .expect("END= present when END_STORED is set");
            self.extract_field((start_key + 4) as i32, consts::SEMICOLON)
        } else {
            String::new()
        }
    }

    /// `toString()` — CHROM, POS, ID, REF/ALT (tab-separated).
    fn to_display(&self) -> String {
        format!(
            "{}{}{}{}{}{}{}",
            self.chrom(),
            consts::TAB,
            self.pos,
            consts::TAB,
            self.id(),
            consts::TAB,
            self.alleles()
        )
    }

    fn append_first_7_fields_and_tab(&self, sb: &mut String) {
        sb.push_str(&self.chrom());
        sb.push(consts::TAB);
        sb.push_str(&self.pos.to_string());
        sb.push(consts::TAB);
        sb.push_str(&self.id());
        sb.push(consts::TAB);
        sb.push_str(&self.alleles());
        sb.push(consts::TAB);
        sb.push_str(&self.qual());
        sb.push(consts::TAB);
        sb.push_str(&self.filter());
        sb.push(consts::TAB);
    }

    /// `appendFirst8Fields(StringBuilder sb)`.
    pub fn append_first_8_fields(&self, sb: &mut String) {
        self.append_first_7_fields_and_tab(sb);
        sb.push_str(&self.info());
    }

    /// `appendFirst8Fields(StringBuilder sb, int an, int[] alleleCounts)` — replaces
    /// INFO/AN and INFO/AC with the supplied counts.
    pub fn append_first_8_fields_with_counts(
        &self,
        sb: &mut String,
        an: i32,
        allele_counts: &[i32],
    ) {
        self.append_first_7_fields_and_tab(sb);
        self.append_info(sb, an, allele_counts);
    }

    fn append_info(&self, sb: &mut String, an: i32, allele_counts: &[i32]) {
        assert!(self.n_alleles() == allele_counts.len() as i32);
        append_counts(sb, an, allele_counts);
        let info = self.info();
        let mut start = 0usize;
        while start < info.len() {
            let mut end = index_of(&info, ';', start);
            if end == -1 {
                end = info.len() as i32;
            }
            if start < end as usize {
                let field = info[start..end as usize].trim();
                if !field.is_empty() && !field.starts_with("AN=") && !field.starts_with("AC=") {
                    sb.push(';');
                    sb.push_str(field);
                }
            }
            start = end as usize + 1;
        }
    }
}

fn append_counts(sb: &mut String, an: i32, allele_counts: &[i32]) {
    sb.push_str("AN=");
    sb.push_str(&an.to_string());
    sb.push_str(";AC=");
    for (j, count) in allele_counts.iter().enumerate().skip(1) {
        if j > 1 {
            sb.push(',');
        }
        sb.push_str(&count.to_string());
    }
}

impl std::fmt::Display for Marker {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.to_display())
    }
}

impl PartialEq for Marker {
    fn eq(&self, other: &Self) -> bool {
        self.chrom_index == other.chrom_index
            && self.pos == other.pos
            && self.alleles() == other.alleles()
            && self.end_value() == other.end_value()
    }
}
impl Eq for Marker {}

impl Hash for Marker {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.chrom_index.hash(state);
        self.pos.hash(state);
        self.alleles().hash(state);
        self.end_value().hash(state);
    }
}

impl Ord for Marker {
    fn cmp(&self, other: &Self) -> Ordering {
        self.chrom_index
            .cmp(&other.chrom_index)
            .then_with(|| self.pos.cmp(&other.pos))
            .then_with(|| self.alleles().cmp(&other.alleles()))
            .then_with(|| self.end_value().cmp(&other.end_value()))
    }
}
impl PartialOrd for Marker {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
