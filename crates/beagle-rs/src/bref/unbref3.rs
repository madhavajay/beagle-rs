//! Port of `bref/UnBref3.java` — the `unbref3` command-line tool: converts a bref3 file back
//! into VCF format, written to standard output.

use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::blbutil::Utilities;
use crate::vcf::{to_vcf_rec, write_meta_lines_gt};

use super::Bref3It;

const PROGRAM: &str = "unbref3.27Feb25.75f.jar";

/// Port of `bref/UnBref3.java`.
pub struct UnBref3;

impl UnBref3 {
    /// `main(String[] args)` — entry point for the `unbref3` tool; writes VCF to stdout.
    pub fn main(args: &[String]) {
        if args.len() > 1 {
            exit(&usage());
        }
        if args.len() == 1 && args[0].eq_ignore_ascii_case("help") {
            exit(&usage());
        }
        let file_name: Option<PathBuf> = if args.is_empty() {
            None
        } else {
            Some(PathBuf::from(&args[0]))
        };
        let mut out = BufWriter::new(std::io::stdout());
        write_vcf(file_name.as_deref(), &mut out);
    }
}

fn exit(msg: &str) -> ! {
    println!("{}", usage());
    Utilities::exit(msg)
}

fn write_vcf<W: Write>(file_name: Option<&Path>, out: &mut W) {
    let mut bref_it = Bref3It::new(file_name);
    if let Some(first) = bref_it.next() {
        let _ = write_meta_lines_gt(&first.samples().ids(), Some(PROGRAM), out);
        let _ = writeln!(out, "{}", to_vcf_rec(first.as_ref()));
    }
    for rec in bref_it.by_ref() {
        let _ = writeln!(out, "{}", to_vcf_rec(rec.as_ref()));
    }
    let _ = out.flush();
}

fn usage() -> String {
    let mut sb = String::with_capacity(500);
    sb.push_str("usage:\n");
    sb.push_str(&format!("  java -jar {PROGRAM} help\n\n"));
    sb.push_str(&format!("  java -jar {PROGRAM} [bref3] > [vcf])\n\n"));
    sb.push_str(&format!(
        "  cat  [bref3]  | java -jar {PROGRAM} > [vcf]\n\n"
    ));
    sb.push_str("where\n");
    sb.push_str("  [bref3]  = the input bref3 file\n");
    sb.push_str("  [vcf]    = the ouput VCF file\n");
    sb
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bref::BrefWriter;
    use crate::bref::CompressBref3Writer;
    use crate::vcf::{
        allele_ref_gt_rec_from_parser, MarkerParser, RefGTRec, VcfHeader, VcfRecGTParser,
        HEADER_PREFIX,
    };
    use std::rc::Rc;

    #[test]
    fn bref_to_vcf_round_trip_lines() {
        // build a small bref3 file on disk, then unbref it to VCF text
        let mut hdr = HEADER_PREFIX.to_string();
        for s in 0..4 {
            hdr.push_str(&format!("\tS{s}"));
        }
        let h = VcfHeader::new_accept_all(
            "src",
            &["##fileformat=VCFv4.2".to_string(), hdr],
            &[true; 4],
        );
        let mp = MarkerParser::new(true, true, true, true);
        let lines = [
            "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\t0|1",
            "chr1\t200\t.\tG\tT\t.\tPASS\t.\tGT\t1|0\t0|0\t0|1\t1|1",
        ];
        let recs: Vec<Rc<dyn RefGTRec>> = lines
            .iter()
            .map(|l| {
                Rc::from(allele_ref_gt_rec_from_parser(&VcfRecGTParser::new(
                    &h, l, &mp,
                )))
            })
            .collect();
        let samples = recs[0].samples().clone();

        let bref_path = std::env::temp_dir().join("beagle_rs_unbref_test.bref3");
        {
            let file = std::fs::File::create(&bref_path).unwrap();
            let mut w = CompressBref3Writer::new(PROGRAM, samples, 7, file);
            for r in &recs {
                w.write(r.clone());
            }
            w.close();
        }

        let mut out: Vec<u8> = Vec::new();
        write_vcf(Some(&bref_path), &mut out);
        let text = String::from_utf8(out).unwrap();
        // meta lines + #CHROM header + 2 records
        assert!(text.contains("##fileformat=VCFv4.2"));
        assert!(text.contains("#CHROM\tPOS\tID"));
        let data_lines: Vec<&str> = text.lines().filter(|l| l.starts_with("chr1\t")).collect();
        assert_eq!(data_lines.len(), 2);
        assert!(data_lines[0].contains("\t100\t"));
        assert!(data_lines[1].contains("\t200\t"));
        let _ = std::fs::remove_file(&bref_path);
    }
}
