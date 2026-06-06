//! Port of `vcf/VcfRecBuilder.java` — builds and writes a VCF 4.2 record line, appending
//! per-sample GT data sequentially.

use crate::blbutil::consts;
use std::io::{self, Write};

use super::{marker_utils, Marker};

/// Port of `vcf/VcfRecBuilder.java`.
pub struct VcfRecBuilder {
    sb: String,
    marker: Marker,
    n_alleles: i32,
}

fn write_fixed_fields(marker: &Marker, sb: &mut String) {
    marker_utils::append_first_7_fields(marker, sb);
    sb.push(consts::TAB);
    sb.push_str(&marker.info());
    sb.push(consts::TAB);
    sb.push_str("GT"); // FORMAT
}

impl VcfRecBuilder {
    /// `new VcfRecBuilder(Marker, int nSamples)`.
    pub fn new(marker: Marker, n_samples: i32) -> Self {
        assert!(n_samples >= 0, "{}", n_samples);
        let mut sb = String::with_capacity(100 + 4 * n_samples.max(0) as usize);
        write_fixed_fields(&marker, &mut sb);
        let n_alleles = marker.n_alleles();
        VcfRecBuilder {
            sb,
            marker,
            n_alleles,
        }
    }

    /// `marker()`.
    pub fn marker(&self) -> &Marker {
        &self.marker
    }

    /// `addSampleData(int a1, int a2)` — diploid phased genotype.
    pub fn add_sample_data_diploid(&mut self, a1: i32, a2: i32) {
        assert!(a1 >= 0 && a1 < self.n_alleles, "{}", a1);
        assert!(a2 >= 0 && a2 < self.n_alleles, "{}", a2);
        self.sb.push(consts::TAB);
        self.sb.push_str(&a1.to_string());
        self.sb.push(consts::PHASED_SEP);
        self.sb.push_str(&a2.to_string());
    }

    /// `addSampleData(int allele)` — haploid phased genotype.
    pub fn add_sample_data_haploid(&mut self, allele: i32) {
        assert!(allele >= 0 && allele < self.n_alleles, "{}", allele);
        self.sb.push(consts::TAB);
        self.sb.push_str(&allele.to_string());
    }

    /// `writeRec(PrintWriter)` — writes the record and a line terminator.
    pub fn write_rec(&self, out: &mut dyn Write) -> io::Result<()> {
        out.write_all(self.sb.as_bytes())?;
        out.write_all(consts::NL.as_bytes())
    }

    /// The record line built so far (without the trailing line terminator).
    pub fn line(&self) -> &str {
        &self.sb
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::MarkerParser;

    #[test]
    fn builds_record_line() {
        let marker = Marker::instance(
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1",
            &MarkerParser::new(true, true, true, true),
        );
        let mut b = VcfRecBuilder::new(marker, 2);
        b.add_sample_data_diploid(0, 1);
        b.add_sample_data_haploid(1);
        let mut out = Vec::new();
        b.write_rec(&mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\t1\n"
        );
    }
}
