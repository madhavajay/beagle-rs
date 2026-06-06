//! Byte-for-byte parity for BGZIP output: the Rust `BgzipOutputStream` must produce
//! exactly the bytes Beagle's `BGZIPOutputStream` produces for the same input
//! (`fixtures/bgzip/sample.bgz`, written by `tools/java/BgzipParityDriver.java`).
//!
//! This is the decisive check that the zlib-backed deflate body matches Java's
//! `Deflater`, which the whole byte-for-byte `.vcf.gz` goal depends on.

use beagle_rs::blbutil::BgzipOutputStream;
use std::io::Write;

fn deterministic_input(n: usize) -> Vec<u8> {
    (0..n).map(|i| ((i * 31 + 7) % 256) as u8).collect()
}

#[test]
fn bgzip_output_matches_java_byte_for_byte() {
    let input = deterministic_input(200_000);
    let mut s = BgzipOutputStream::new(Vec::new(), true);
    s.write_all(&input).unwrap();
    let out = s.close().unwrap();
    let expected: &[u8] = include_bytes!("../../../fixtures/bgzip/sample.bgz");
    assert_eq!(
        out.len(),
        expected.len(),
        "compressed length differs from Java"
    );
    assert_eq!(
        out, expected,
        "BGZIP output diverged from the Java reference"
    );
}
