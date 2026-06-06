//! Port of `imp/ImputedRecBuilder.java` — accumulates per-sample imputed allele probabilities
//! and prints a VCF 4.3 record (GT:DS[:AP1:AP2][:GP]) with a DR2/AF INFO field.
//!
//! Byte-exact formatting: the Java `DecimalFormat` tables `DS_VALS` ("#.##") and `R2_VALS`
//! ("0.00") reduce to integer-cents formatting (verified against Java 26 output); the AF
//! subfield uses Java's `DecimalFormat("0.0000")`, reproduced by Rust's `{:.4}` (both round
//! half-to-even on the exact `f64`).

use std::io::Write;

use crate::blbutil::consts;
use crate::vcf::Marker;

const DEFAULT_HOM_REF_LEN: i32 = 5;

/// `DS_VALS[index]` — Java `new DecimalFormat("#.##").format(index/100.0)` for `index` in
/// `0..=200`: integer-cents with trailing zeros dropped.
fn ds_val(index: i32) -> String {
    let whole = index / 100;
    let frac = index % 100;
    if frac == 0 {
        whole.to_string()
    } else if frac % 10 == 0 {
        format!("{}.{}", whole, frac / 10)
    } else {
        format!("{}.{:02}", whole, frac)
    }
}

/// `R2_VALS[index]` — Java `new DecimalFormat("0.00").format(index/100.0)` for `index` in
/// `0..=100`: always two fractional digits.
fn r2_val(index: i32) -> String {
    format!("{}.{:02}", index / 100, index % 100)
}

/// `(int) Math.rint(100 * x)` — `x` is `float`; `100*x` is computed in `f32`, widened to
/// `f64` for round-half-even, then truncated.
fn rint100(x: f32) -> i32 {
    ((100.0f32 * x) as f64).round_ties_even() as i32
}

fn scale(fa: &mut [f32]) {
    let mut sum = 0f32;
    for &f in fa.iter() {
        sum += f;
    }
    for f in fa.iter_mut() {
        *f /= sum;
    }
}

fn max_index(fa: &[f32]) -> i32 {
    let mut max_index = 0usize;
    for j in 1..fa.len() {
        if fa[j] > fa[max_index] {
            max_index = j;
        }
    }
    max_index as i32
}

fn default_hom_ref_fields() -> Vec<String> {
    let mut sa = vec![String::new(); DEFAULT_HOM_REF_LEN as usize];
    sa[1] = format!("{}0|0", consts::TAB);
    sa[2] = format!("{}0|0:0", consts::TAB);
    for j in 3..sa.len() {
        sa[j] = format!("{},0", sa[j - 1]);
    }
    sa
}

fn hom_ref_fields(ap: bool, gp: bool) -> Vec<String> {
    let default = default_hom_ref_fields();
    let mut hom_ref_field = vec![String::new(); default.len()];
    for n_al in 1..hom_ref_field.len() {
        let mut sb = default[n_al].clone();
        if ap {
            for a in 1..n_al {
                sb.push(if a == 1 { consts::COLON } else { consts::COMMA });
                sb.push_str(&ds_val(0));
            }
            for a in 1..n_al {
                sb.push(if a == 1 { consts::COLON } else { consts::COMMA });
                sb.push_str(&ds_val(0));
            }
        }
        if gp {
            sb.push(consts::COLON);
            sb.push_str(&ds_val(100));
            for i2 in 1..n_al {
                for _i1 in 0..=i2 {
                    sb.push(consts::COMMA);
                    sb.push_str(&ds_val(0));
                }
            }
        }
        hom_ref_field[n_al] = sb;
    }
    hom_ref_field
}

/// Port of `imp/ImputedRecBuilder.java`.
pub struct ImputedRecBuilder {
    marker: Marker,
    n_alleles: i32,
    n_input_targ_haps: i32,
    ap: bool,
    gp: bool,
    sum_al_probs: Vec<f32>,
    sum_al_probs2: Vec<f32>,
    hom_ref_field: Vec<String>,
    sample_data: String,
    hap_cnt: i32,
}

impl ImputedRecBuilder {
    /// `new ImputedRecBuilder(Marker, int nInputTargHaps, boolean ap, boolean gp)`.
    pub fn new(marker: Marker, n_input_targ_haps: i32, ap: bool, gp: bool) -> Self {
        assert!(n_input_targ_haps >= 1, "{n_input_targ_haps}");
        let n_alleles = marker.n_alleles();
        let hom_ref_field = if ap || gp {
            hom_ref_fields(ap, gp)
        } else {
            default_hom_ref_fields()
        };
        ImputedRecBuilder {
            marker,
            n_alleles,
            n_input_targ_haps,
            ap,
            gp,
            sum_al_probs: vec![0.0; n_alleles as usize],
            sum_al_probs2: vec![0.0; n_alleles as usize],
            hom_ref_field,
            sample_data: String::with_capacity(200 + n_input_targ_haps as usize * 5),
            hap_cnt: 0,
        }
    }

    /// `marker()`.
    pub fn marker(&self) -> &Marker {
        &self.marker
    }
    /// `nInputTargHaps()`.
    pub fn n_input_targ_haps(&self) -> i32 {
        self.n_input_targ_haps
    }
    /// `hapCnt()`.
    pub fn hap_cnt(&self) -> i32 {
        self.hap_cnt
    }

    /// `addSampleData(float[] a1, float[] a2)` — diploid sample.
    // `a` indexes the allele probabilities and drives the `a==1` field-vs-subfield separator.
    #[allow(clippy::needless_range_loop)]
    pub fn add_sample_data(&mut self, a1: &mut [f32], a2: &mut [f32]) {
        self.hap_cnt += 2;
        if a1[0] == 1.0 && a2[0] == 1.0 && (a1.len() as i32) < DEFAULT_HOM_REF_LEN {
            self.sample_data.push_str(&self.hom_ref_field[a1.len()]);
            return;
        }
        scale(a1);
        scale(a2);
        self.sample_data.push(consts::TAB);
        self.sample_data.push_str(&max_index(a1).to_string());
        self.sample_data.push(consts::PHASED_SEP);
        self.sample_data.push_str(&max_index(a2).to_string());
        let n = self.n_alleles as usize;
        for a in 1..n {
            let dose = a1[a] + a2[a];
            let dose2 = a1[a] * a1[a] + a2[a] * a2[a];
            self.sum_al_probs[a] += dose;
            self.sum_al_probs2[a] += dose2;
            self.sample_data
                .push(if a == 1 { consts::COLON } else { consts::COMMA });
            self.sample_data.push_str(&ds_val(rint100(dose)));
        }
        if self.ap {
            for a in 1..n {
                self.sample_data
                    .push(if a == 1 { consts::COLON } else { consts::COMMA });
                self.sample_data.push_str(&ds_val(rint100(a1[a])));
            }
            for a in 1..n {
                self.sample_data
                    .push(if a == 1 { consts::COLON } else { consts::COMMA });
                self.sample_data.push_str(&ds_val(rint100(a2[a])));
            }
        }
        if self.gp {
            for i2 in 0..n {
                for i1 in 0..=i2 {
                    let mut prob = a1[i1] * a2[i2];
                    if i1 != i2 {
                        prob += a1[i2] * a2[i1];
                    }
                    self.sample_data.push(if i2 == 0 {
                        consts::COLON
                    } else {
                        consts::COMMA
                    });
                    self.sample_data.push_str(&ds_val(rint100(prob)));
                }
            }
        }
    }

    /// `addSampleData(float[] a1)` — haploid sample.
    #[allow(clippy::needless_range_loop)]
    pub fn add_sample_data_haploid(&mut self, a1: &mut [f32]) {
        self.hap_cnt += 1;
        scale(a1);
        self.sample_data.push(consts::TAB);
        self.sample_data.push_str(&max_index(a1).to_string());
        let n = self.n_alleles as usize;
        for a in 1..n {
            let dose = a1[a];
            let dose2 = a1[a] * a1[a];
            self.sum_al_probs[a] += dose;
            self.sum_al_probs2[a] += dose2;
            self.sample_data
                .push(if a == 1 { consts::COLON } else { consts::COMMA });
            self.sample_data.push_str(&ds_val(rint100(dose)));
        }
        if self.ap {
            for a in 1..n {
                self.sample_data
                    .push(if a == 1 { consts::COLON } else { consts::COMMA });
                self.sample_data.push_str(&ds_val(rint100(a1[a])));
            }
        }
    }

    fn r2(&self, allele: usize) -> f32 {
        let sum = self.sum_al_probs[allele];
        if sum == 0.0 {
            0.0
        } else {
            let sum2 = self.sum_al_probs2[allele];
            let mean_term = sum * sum / self.n_input_targ_haps as f32;
            let num = sum2 - mean_term;
            let den = sum - mean_term;
            if num <= 0.0 {
                0.0
            } else {
                num / den
            }
        }
    }

    /// `printRec(PrintWriter out, boolean isImputed)`.
    pub fn print_rec(&self, out: &mut dyn Write, is_imputed: bool) {
        assert!(self.hap_cnt == self.n_input_targ_haps, "inconsistent data");
        let m = &self.marker;
        let _ = write!(
            out,
            "{}{tab}{}{tab}{}{tab}{}",
            m.chrom(),
            m.pos(),
            m.id(),
            m.alleles(),
            tab = consts::TAB
        );
        let _ = write!(
            out,
            "{tab}{miss}{tab}PASS{tab}",
            tab = consts::TAB,
            miss = consts::MISSING_DATA_CHAR
        );
        self.print_info_field(out, is_imputed);
        let _ = write!(out, "{}GT:DS", consts::TAB);
        if self.ap {
            let _ = write!(out, ":AP1:AP2");
        }
        if self.gp {
            let _ = write!(out, ":GP");
        }
        let _ = writeln!(out, "{}", self.sample_data);
    }

    fn print_info_field(&self, out: &mut dyn Write, is_imputed: bool) {
        if self.n_alleles == 1 {
            if is_imputed {
                let _ = write!(out, "IMP");
            }
            return;
        }
        for a in 1..self.n_alleles as usize {
            let _ = write!(out, "{}", if a == 1 { "DR2=" } else { "," });
            let _ = write!(out, "{}", r2_val(rint100(self.r2(a))));
        }
        for a in 1..self.n_alleles as usize {
            let _ = write!(out, "{}", if a == 1 { ";AF=" } else { "," });
            let af = self.sum_al_probs[a] / self.n_input_targ_haps as f32;
            let _ = write!(out, "{:.4}", af as f64);
        }
        if let Some(end_subfield) = extract_end(&self.marker) {
            let _ = write!(out, ";{}", end_subfield);
        }
        if is_imputed {
            let _ = write!(out, ";IMP");
        }
    }
}

fn extract_end(marker: &Marker) -> Option<String> {
    let info = marker.info();
    let start = info.find("END=")?;
    let end_index = info[start + 4..]
        .find(consts::SEMICOLON)
        .map(|p| start + 4 + p);
    let e = end_index.unwrap_or(info.len());
    Some(info[start..e].to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::MarkerParser;

    fn mk(rec: &str) -> Marker {
        Marker::instance(rec, &MarkerParser::new(true, true, true, true))
    }

    #[test]
    fn ds_and_r2_tables_match_java() {
        // verified against Java DecimalFormat("#.##")/("0.00")
        assert_eq!(ds_val(0), "0");
        assert_eq!(ds_val(5), "0.05");
        assert_eq!(ds_val(10), "0.1");
        assert_eq!(ds_val(50), "0.5");
        assert_eq!(ds_val(100), "1");
        assert_eq!(ds_val(150), "1.5");
        assert_eq!(ds_val(200), "2");
        assert_eq!(ds_val(7), "0.07");
        assert_eq!(r2_val(0), "0.00");
        assert_eq!(r2_val(50), "0.50");
        assert_eq!(r2_val(100), "1.00");
        assert_eq!(r2_val(5), "0.05");
    }

    #[test]
    fn biallelic_record_round_trip() {
        let marker = mk("chr1\t100\trs1\tA\tC\t.\tPASS\t.\tGT");
        let mut b = ImputedRecBuilder::new(marker, 4, false, false);
        // two diploid samples
        b.add_sample_data(&mut [1.0, 0.0], &mut [1.0, 0.0]); // hom ref -> homRefField[2]
        b.add_sample_data(&mut [0.0, 1.0], &mut [0.5, 0.5]);
        assert_eq!(b.hap_cnt(), 4);
        let mut out: Vec<u8> = Vec::new();
        b.print_rec(&mut out, true);
        let text = String::from_utf8(out).unwrap();
        assert!(text.starts_with("chr1\t100\trs1\tA\tC\t.\tPASS\t"));
        assert!(text.contains("DR2="));
        assert!(text.contains(";AF="));
        assert!(text.contains(";IMP"));
        assert!(text.contains("GT:DS\t"));
        // first sample is hom-ref "0|0:0"
        assert!(text.contains("\t0|0:0"));
        assert!(text.ends_with('\n'));
    }

    #[test]
    fn monomorphic_record_info_is_imp_only() {
        let marker = mk("chr1\t100\t.\tA\t.\t.\tPASS\t.\tGT");
        let mut b = ImputedRecBuilder::new(marker, 2, false, false);
        b.add_sample_data(&mut [1.0], &mut [1.0]);
        let mut out: Vec<u8> = Vec::new();
        b.print_rec(&mut out, true);
        let text = String::from_utf8(out).unwrap();
        // INFO is just IMP for a monomorphic marker
        assert!(text.contains("\tIMP\tGT:DS\t"));
    }
}
