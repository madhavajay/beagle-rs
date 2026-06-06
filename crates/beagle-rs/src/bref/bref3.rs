//! Port of `bref/Bref3.java` — the `bref3` command-line tool: converts a VCF file (phased,
//! non-missing genotypes) into bref3 format, written to standard output.

use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::rc::Rc;

use crate::blbutil::{InputIt, SampleFileIt, Utilities};
use crate::vcf::RefIt;

use super::{AsIsBref3Writer, BrefWriter, CompressBref3Writer};

const PROGRAM: &str = "bref3.27Feb25.75f.jar";
const CHARACTER_MAX_VALUE: i32 = 65535;

/// Port of `bref/Bref3.java`.
pub struct Bref3;

impl Bref3 {
    /// `main(String[] args)` — entry point for the `bref3` tool; writes bref3 to stdout.
    pub fn main(args: &[String]) {
        if args.len() > 2 {
            println!("{}", usage());
            std::process::exit(0);
        }
        if args.len() == 1 && args[0].eq_ignore_ascii_case("help") {
            println!("{}", usage());
            std::process::exit(0);
        }
        let use_std_in = use_std_in(args);
        let mut max_n_seq = -1;
        if args.len() == 2 || (use_std_in && args.len() == 1) {
            max_n_seq = parse_max_n_seq(&args[args.len() - 1]);
        }
        let input_file: Option<PathBuf> = if use_std_in {
            None
        } else {
            Some(PathBuf::from(&args[0]))
        };
        write_bref(input_file.as_deref(), max_n_seq);
    }
}

fn use_std_in(sa: &[String]) -> bool {
    match sa.len() {
        0 => true,
        1 => !(sa[0].ends_with(".vcf") || sa[0].ends_with(".vcf.gz")),
        _ => false,
    }
}

fn parse_max_n_seq(arg: &str) -> i32 {
    match arg.parse::<i32>() {
        Ok(n) => {
            if !(1..=CHARACTER_MAX_VALUE).contains(&n) {
                exit(&format!("Error: invalid <nSeq> {arg}"));
            }
            n
        }
        Err(_) => exit(&format!("Error: <nSeq> is not a parsable integer: {arg}")),
    }
}

fn exit(msg: &str) -> ! {
    println!("{}", usage());
    Utilities::exit(msg)
}

fn ref_it(file_name: Option<&Path>) -> RefIt<InputIt> {
    match file_name {
        None => RefIt::create(InputIt::from_reader(
            Box::new(BufReader::new(std::io::stdin())),
            None,
        )),
        Some(file) => {
            let n_cores = std::thread::available_parallelism()
                .map(|n| n.get() as i32)
                .unwrap_or(1);
            let n_buffered_blocks = n_cores << 2;
            RefIt::create(InputIt::from_bgzip_file(file, n_buffered_blocks))
        }
    }
}

fn write_bref(file_name: Option<&Path>, max_n_seq: i32) {
    let mut it = ref_it(file_name);
    let samples = it.samples().clone();
    let out: Box<dyn Write> = Box::new(BufWriter::new(std::io::stdout()));
    let mut bref_out: Box<dyn BrefWriter> = if max_n_seq < 0 {
        Box::new(AsIsBref3Writer::new(PROGRAM, samples, out))
    } else {
        Box::new(CompressBref3Writer::new(PROGRAM, samples, max_n_seq, out))
    };
    for rec in it.by_ref() {
        bref_out.write(Rc::from(rec));
    }
    bref_out.close();
}

fn usage() -> String {
    let mut sb = String::with_capacity(500);
    sb.push_str("usage:\n");
    sb.push_str(&format!("  java -jar {PROGRAM} help\n\n"));
    sb.push_str(&format!(
        "  java -jar {PROGRAM} [vcf] <nseq>  > [bref3]\n\n"
    ));
    sb.push_str(&format!(
        "  cat   [vcf]   | java -jar {PROGRAM} <nseq>  > [bref3]\n\n"
    ));
    sb.push_str("where\n");
    sb.push_str("  [bref3]  = the output bref3 file\n");
    sb.push_str("  [vcf]    = A VCF file with phased, non-missing genotype data.  If the\n");
    sb.push_str("             file is gzip-compressed, its filename must end in \".gz\"\n");
    sb.push_str("             and \"cat\" must be replaced with \"zcat\"\n");
    sb.push_str("  <nseq>   = optional argument for maximum number of unique sequences\n");
    sb.push_str("             in a bref3 block. If there are N reference samples,\n");
    sb.push_str("             the default value is: <max-seq>=2^(2*log10(N) + 1)\n");
    sb
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn use_std_in_and_usage() {
        assert!(use_std_in(&[]));
        assert!(use_std_in(&["5".to_string()])); // not a vcf name -> nSeq via stdin
        assert!(!use_std_in(&["x.vcf".to_string()]));
        assert!(!use_std_in(&["x.vcf.gz".to_string()]));
        assert!(!use_std_in(&["a".to_string(), "b".to_string()]));
        assert!(usage().contains("bref3"));
    }

    #[test]
    fn parses_max_n_seq() {
        assert_eq!(parse_max_n_seq("100"), 100);
        assert_eq!(parse_max_n_seq("1"), 1);
    }
}
