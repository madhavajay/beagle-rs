//! Cross-language parity for `vcf::Marker` (+ `MarkerUtils`, `MarkerParser`). The Rust
//! port must reproduce `fixtures/marker/parity.txt` (from
//! `tools/java/MarkerParityDriver.java`) — parsing real VCF records and matching every
//! accessor, for both store-all and store-none parser configs.

use beagle_rs::vcf::{Marker, MarkerParser};

const RECORDS: &[&str] = &[
    "chr1\t100\trs1\tA\tC\t30\tPASS\t.\tGT\t0|1",
    "chr1\t200\t.\tA\tG,T\t.\t.\tAC=2;AN=4\tGT\t1|2",
    "chr1\t300\tidDEL\tAC\tA\t.\tq10\tEND=305;SVTYPE=DEL\tGT\t0|0",
    "chr1\t400\t.\tA\t.\t.\t.\t.\tGT\t0|0",
    "chr2\t500\trs5\tACGT\tA,ATTT\t99\tPASS\tDP=10\tGT\t0|1",
    "chr2\t600\t.\tG\tA,C,T\t.\tPASS\t.\tGT\t0|3",
];

fn clean(s: &str) -> String {
    s.replace('\t', "/")
}

fn emit(out: &mut String, mode: &str, parser: &MarkerParser) {
    for (i, rec) in RECORDS.iter().enumerate() {
        let m = Marker::instance(rec, parser);
        out.push_str(&format!(
            "mode={} rec={} chrom={} pos={} id={} alleles={} nAlleles={} nRefBases={} qual={} \
             filter={} info={} end={} hasId={} hasEnd={} bits={} str={}\n",
            mode,
            i,
            m.chrom(),
            m.pos(),
            m.id(),
            clean(&m.alleles()),
            m.n_alleles(),
            m.n_ref_bases(),
            m.qual(),
            m.filter(),
            m.info(),
            m.end_value(),
            m.has_id_data(),
            m.has_end_value(),
            m.bits_per_allele(),
            clean(&m.to_string()),
        ));
    }
}

fn transcript() -> String {
    let mut out = String::new();
    out.push_str("# marker parity (Beagle 5.5 27Feb25.75f)\n");
    emit(&mut out, "full", &MarkerParser::new(true, true, true, true));
    emit(
        &mut out,
        "none",
        &MarkerParser::new(false, false, false, false),
    );
    out
}

#[test]
fn marker_parity_matches_java_reference() {
    let expected = include_str!("../../../fixtures/marker/parity.txt");
    assert_eq!(
        transcript(),
        expected,
        "Rust Marker accessors diverged from the Java reference"
    );
}
