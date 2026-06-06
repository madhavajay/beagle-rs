//! Port of `vcf/MarkerUtils.java` — static helpers for `Marker`, including the sorted
//! SNV-permutation table used to compactly encode single-nucleotide/monomorphic alleles.

use crate::beagleutil::ChromIds;
use crate::blbutil::{consts, StringUtil, Utilities};
use std::sync::OnceLock;

use super::Marker;

/// `Marker.SNV_PERMS` — sorted REF/ALT strings for SNV and monomorphic variants. Each
/// non-monomorphic entry is the maximal permutation `R\tA,B,C`; `Marker.alleles()`
/// takes a prefix for fewer alleles.
pub(crate) fn snv_perms() -> &'static Vec<String> {
    static PERMS: OnceLock<Vec<String>> = OnceLock::new();
    PERMS.get_or_init(|| {
        let mut bases = ['*', 'A', 'C', 'G', 'T'];
        bases.sort_unstable();
        let mut perms: Vec<String> = Vec::with_capacity(120);
        permute(&[], &bases, &mut perms);
        perms.push("A\t.".to_string());
        perms.push("C\t.".to_string());
        perms.push("G\t.".to_string());
        perms.push("T\t.".to_string());
        perms.sort();
        perms
    })
}

fn permute(start: &[char], end: &[char], perms: &mut Vec<String>) {
    if end.is_empty() && start[0] != '*' {
        let mut s = String::new();
        s.push(start[0]);
        for (j, &c) in start.iter().enumerate().skip(1) {
            s.push(if j == 1 { consts::TAB } else { consts::COMMA });
            s.push(c);
        }
        perms.push(s);
    } else {
        for j in 0..end.len() {
            let mut new_start = start.to_vec();
            new_start.push(end[j]);
            let mut new_end = Vec::with_capacity(end.len() - 1);
            new_end.extend_from_slice(&end[0..j]);
            new_end.extend_from_slice(&end[j + 1..]);
            permute(&new_start, &new_end, perms);
        }
    }
}

/// `nGenotypes(int nAlleles)` = `nAlleles*(nAlleles+1)/2`.
pub fn n_genotypes(n_alleles: i32) -> i32 {
    (n_alleles * (n_alleles + 1)) >> 1
}

/// `truncate(String s, int maxLength)`.
pub(crate) fn truncate(s: &str, max_length: usize) -> &str {
    &s[0..s.len().min(max_length)]
}

/// `first8TabIndices(String vcfRec)` — byte indices of the first 8 tabs.
pub(crate) fn first_8_tab_indices(vcf_rec: &str) -> Vec<i32> {
    const N_TABS: usize = 8;
    let mut indices = Vec::with_capacity(N_TABS);
    for (i, b) in vcf_rec.bytes().enumerate() {
        if b == b'\t' {
            indices.push(i as i32);
            if indices.len() == N_TABS {
                break;
            }
        }
    }
    if indices.len() < N_TABS {
        panic!(
            "VCF record does not contain {N_TABS} tabs:{}",
            truncate(vcf_rec, 800)
        );
    }
    indices
}

/// `VcfRecGTParser.ninthTabPos(vcfRec)` — byte index of the 9th tab (end of FORMAT).
pub(crate) fn ninth_tab_pos(vcf_rec: &str) -> i32 {
    let bytes = vcf_rec.as_bytes();
    let mut pos: i32 = -1;
    for _ in 0..9 {
        match bytes[(pos + 1) as usize..].iter().position(|&b| b == b'\t') {
            Some(off) => pos = (pos + 1) + off as i32,
            None => panic!("VCF record format error: {vcf_rec}"),
        }
    }
    pos
}

/// `chromIndex(String vcfRec, String chrom)` — validates and indexes the CHROM field.
pub(crate) fn chrom_index(vcf_rec: &str, chrom: &str) -> i16 {
    if chrom.is_empty() || chrom == "." {
        panic!("ERROR: missing chromosome: {}", truncate(vcf_rec, 80));
    }
    for c in chrom.chars() {
        if c.is_whitespace() {
            Utilities::exit(&format!(
                "ERROR: CHROM field contains whitespace: {}",
                truncate(vcf_rec, 80)
            ));
        }
    }
    let chr_index = ChromIds::instance().get_index(chrom);
    assert!(chr_index < i16::MAX as i32, "{}", chr_index);
    chr_index as i16
}

/// `appendFirst7Fields(Marker, StringBuilder)` — CHROM..FILTER, no trailing tab.
pub fn append_first_7_fields(marker: &Marker, sb: &mut String) {
    sb.push_str(&marker.chrom());
    sb.push(consts::TAB);
    sb.push_str(&marker.pos().to_string());
    sb.push(consts::TAB);
    sb.push_str(&marker.id());
    sb.push(consts::TAB);
    sb.push_str(&marker.alleles());
    sb.push(consts::TAB);
    sb.push_str(&marker.qual());
    sb.push(consts::TAB);
    sb.push_str(&marker.filter());
}

/// `coordinate(Marker)` — `chrom:pos`.
pub fn coordinate(marker: &Marker) -> String {
    format!("{}{}{}", marker.chrom(), consts::COLON, marker.pos())
}

/// `coordinateAndAlleles(Marker)` — `chrom:pos:REF:ALT`.
pub fn coordinate_and_alleles(marker: &Marker) -> String {
    format!(
        "{}{}{}{}{}",
        marker.chrom(),
        consts::COLON,
        marker.pos(),
        consts::COLON,
        marker
            .alleles()
            .replace(consts::TAB, &consts::COLON.to_string())
    )
}

/// `ids(Marker)` — the `;`-separated identifiers in the VCF ID field (empty if none).
pub fn ids(marker: &Marker) -> Vec<String> {
    if marker.has_id_data() {
        let id = marker.id();
        StringUtil::get_fields(&id, consts::SEMICOLON as u8)
            .into_iter()
            .map(|s| s.to_string())
            .collect()
    } else {
        Vec::new()
    }
}

/// `alleles(Marker)` — the alleles (REF then ALTs) in VCF order.
pub fn alleles(marker: &Marker) -> Vec<String> {
    let ref_alt = marker.alleles();
    let n_alleles = marker.n_alleles();
    let mut out = Vec::with_capacity(n_alleles as usize);
    let tab = ref_alt.find('\t').unwrap();
    out.push(ref_alt[0..tab].to_string());
    let mut start = tab + 1;
    for _ in 1..n_alleles {
        match ref_alt[start..].find(',') {
            Some(off) => {
                let end = start + off;
                out.push(ref_alt[start..end].to_string());
                start = end + 1;
            }
            None => {
                out.push(ref_alt[start..].to_string());
                start = ref_alt.len() + 1;
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snv_perms_table() {
        let perms = snv_perms();
        assert_eq!(perms.len(), 100); // matches the Java comment
        assert!(perms.is_sorted());
        // every monomorphic form present
        for m in ["A\t.", "C\t.", "G\t.", "T\t."] {
            assert!(perms.iter().any(|p| p == m));
        }
        // a known biallelic prefix exists as a maximal perm
        assert!(perms.iter().any(|p| p.starts_with("A\tC")));
    }

    #[test]
    fn n_genotypes_and_truncate() {
        assert_eq!(n_genotypes(2), 3);
        assert_eq!(n_genotypes(3), 6);
        assert_eq!(truncate("hello", 3), "hel");
        assert_eq!(truncate("hi", 10), "hi");
    }

    #[test]
    fn tab_indices() {
        let rec = "chr1\t100\trs1\tA\tC\t.\tPASS\t.\tGT\t0|1";
        assert_eq!(first_8_tab_indices(rec).len(), 8);
    }

    #[test]
    #[should_panic]
    fn tab_indices_too_few() {
        let _ = first_8_tab_indices("chr1\t100\tonly\ttwo\ttabs");
    }
}
