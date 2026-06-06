//! Port of `vcf/MarkerParser.java` — parses and filters a VCF record's ID/REF/ALT/
//! QUAL/FILTER/INFO subfields into the packed `Marker` `field_info` + `fields` string.

use crate::blbutil::{consts, StringUtil};

use super::marker::{
    ALLELES_STORED, END_STORED, FILTER_STORED, ID_STORED, INFO_STORED, QUAL_STORED,
    STORED_N_ALLELES_MASK,
};
use super::marker_utils;

/// `String.indexOf(ch, from)` — byte index or -1.
fn index_of(s: &str, ch: char, from: usize) -> i32 {
    s[from..].find(ch).map_or(-1, |i| (from + i) as i32)
}

/// Port of `vcf/MarkerParser.java`.
#[derive(Clone, Copy, Debug)]
pub struct MarkerParser {
    store_id: bool,
    store_qual: bool,
    store_filter: bool,
    store_info: bool,
}

impl MarkerParser {
    /// `new MarkerParser(storeId, storeQual, storeFilter, storeInfo)`.
    pub fn new(store_id: bool, store_qual: bool, store_filter: bool, store_info: bool) -> Self {
        MarkerParser {
            store_id,
            store_qual,
            store_filter,
            store_info,
        }
    }

    /// `storeMarkerFields(String rec, StringBuilder sb, IntList tabs)`.
    pub(crate) fn store_marker_fields(&self, rec: &str, sb: &mut String, tabs: &[i32]) -> u16 {
        let mut info: u16 = 0;
        info = store_field(
            rec,
            tabs[1] + 1,
            tabs[2],
            info,
            sb,
            self.store_id,
            ID_STORED,
        );
        info = store_alleles(rec, tabs[2] + 1, tabs[4], info, sb);
        info = store_field(
            rec,
            tabs[4] + 1,
            tabs[5],
            info,
            sb,
            self.store_qual,
            QUAL_STORED,
        );
        info = store_field(
            rec,
            tabs[5] + 1,
            tabs[6],
            info,
            sb,
            self.store_filter,
            FILTER_STORED,
        );
        info = self.store_info(rec, tabs[6] + 1, tabs[7], info, sb);
        info
    }

    fn store_info(
        &self,
        vcf_rec: &str,
        start: i32,
        end: i32,
        mut flags: u16,
        sb: &mut String,
    ) -> u16 {
        assert!(
            start >= 0 && start <= end && end <= vcf_rec.len() as i32,
            "start={start} end={end}"
        );
        let (s, e) = (start as usize, end as usize);
        if (end - start) != 1 || vcf_rec.as_bytes()[s] != b'.' {
            let end_subfield = end_subfield(vcf_rec, s, e);
            if self.store_info || end_subfield.is_some() {
                if !sb.is_empty() {
                    sb.push(consts::TAB);
                }
                if self.store_info {
                    sb.push_str(&vcf_rec[s..e]);
                    flags |= INFO_STORED;
                    if end_subfield.is_some() {
                        flags |= END_STORED;
                    }
                } else {
                    sb.push_str(&end_subfield.expect("END subfield present"));
                    flags |= END_STORED;
                }
            }
        }
        flags
    }

    /// `storeId()`.
    pub fn store_id(&self) -> bool {
        self.store_id
    }
    /// `storeQual()`.
    pub fn store_qual(&self) -> bool {
        self.store_qual
    }
    /// `storeFilter()`.
    pub fn store_filter(&self) -> bool {
        self.store_filter
    }
    /// `storeInfo()`.
    pub fn store_info_flag(&self) -> bool {
        self.store_info
    }
}

#[allow(clippy::too_many_arguments)]
fn store_field(
    vcf_rec: &str,
    start: i32,
    end: i32,
    mut flags: u16,
    sb: &mut String,
    is_stored: bool,
    flag: u16,
) -> u16 {
    assert!(
        start >= 0 && start <= end && end <= vcf_rec.len() as i32,
        "start={start} end={end}"
    );
    let (s, e) = (start as usize, end as usize);
    if is_stored && ((end - start) != 1 || vcf_rec.as_bytes()[s] != b'.') {
        if !sb.is_empty() {
            sb.push(consts::TAB);
        }
        sb.push_str(&vcf_rec[s..e]);
        flags |= flag;
    }
    flags
}

fn store_alleles(vcf_rec: &str, start: i32, end: i32, mut field_info: u16, sb: &mut String) -> u16 {
    let ref_alt = &vcf_rec[start as usize..end as usize];
    let snv_index = snv_index(ref_alt);
    if snv_index >= 0 {
        let n_alleles = if ref_alt.ends_with("\t.") {
            1
        } else {
            (end - start + 1) >> 1
        };
        field_info |= (snv_index as u16) << 3;
        field_info |= n_alleles as u16;
    } else {
        let tab_index = index_of(ref_alt, consts::TAB, 0);
        if (tab_index + 1) as usize == ref_alt.len() {
            panic!(
                "ERROR: missing ALT field: {}",
                marker_utils::truncate(vcf_rec, 80)
            );
        }
        let n_alleles = n_alleles_of(ref_alt, tab_index);
        if n_alleles > STORED_N_ALLELES_MASK as i32 {
            panic!(
                "{} alleles: {}",
                n_alleles,
                &vcf_rec[0..80.min(vcf_rec.len())]
            );
        }
        if !sb.is_empty() {
            sb.push(consts::TAB);
        }
        sb.push_str(ref_alt);
        field_info |= n_alleles as u16;
        field_info |= ALLELES_STORED;
    }
    field_info
}

fn snv_index(ref_and_alt: &str) -> i32 {
    let perms = marker_utils::snv_perms();
    let index = match perms.binary_search_by(|p| p.as_str().cmp(ref_and_alt)) {
        Ok(i) | Err(i) => i,
    };
    if index == perms.len() {
        -1
    } else if perms[index].starts_with(ref_and_alt) {
        index as i32
    } else {
        -1
    }
}

fn n_alleles_of(ref_alt_alleles: &str, tab_index: i32) -> i32 {
    let start_allele = tab_index + 1;
    let bytes = ref_alt_alleles.as_bytes();
    if start_allele as usize == ref_alt_alleles.len() - 1
        && bytes[start_allele as usize] == consts::MISSING_DATA_CHAR as u8
    {
        1
    } else {
        let mut n_alleles = 2;
        let mut sa = index_of(ref_alt_alleles, consts::COMMA, start_allele as usize) + 1;
        while sa > 0 {
            n_alleles += 1;
            sa = index_of(ref_alt_alleles, consts::COMMA, sa as usize) + 1;
        }
        n_alleles
    }
}

fn end_subfield(vcf_rec: &str, info_start: usize, info_end: usize) -> Option<String> {
    let info_field = &vcf_rec[info_start..info_end];
    for field in StringUtil::get_fields(info_field, consts::SEMICOLON as u8) {
        if field.starts_with("END=") {
            return Some(field.to_string());
        }
    }
    None
}

impl std::fmt::Display for MarkerParser {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "ID={} QUAL={} FILTER={} INFO={}",
            self.store_id, self.store_qual, self.store_filter, self.store_info
        )
    }
}
