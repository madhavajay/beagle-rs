//! End-to-end byte-for-byte parity against the Java reference (`beagle.27Feb25.75f.jar`)
//! on the committed official fixtures (`fixtures/official/`).
//!
//! These run the actual CLI binaries on the official 1356-marker × 191-sample slice and assert
//! the produced `*.vcf.gz` / `*.bref3` is byte-identical to the golden output the Java reference
//! generated with `nthreads=1 seed=99999` (see `fixtures/official/MANIFEST.md`).
//!
//! The golden `.vcf.gz` embeds a wall-clock `##filedate` (the documented normalization field). To
//! make the byte comparison date-independent *and still prove BGZIP-byte parity*, the runs set
//! `SOURCE_DATE_EPOCH` to an instant on the golden's date (2026-06-06 UTC) so `##filedate` matches
//! the committed golden exactly — then every byte, including the deflate stream, must be identical.
//!
//! Marked `#[ignore]` because the `gt=` phasing run takes ~16 s in a debug build (≈2 s in
//! release). CI runs them explicitly in release:
//!   `cargo test -p beagle-rs-cli --release --test parity_official -- --ignored`
//! or via `scripts/parity-check.sh`. They skip gracefully if the fixtures are absent.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// An instant on 2026-06-06 UTC — the date baked into the committed golden `##filedate`. Passed
/// via `SOURCE_DATE_EPOCH` so the wall-clock field matches and the comparison is byte-exact.
const GOLDEN_SOURCE_DATE_EPOCH: &str = "1780747200"; // 2026-06-06 12:00:00 UTC

fn fixtures_dir() -> Option<PathBuf> {
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/official");
    (d.join("golden").is_dir() && d.join("inputs").is_dir()).then_some(d)
}

/// Copy `fixtures/official/inputs/*` into a fresh temp dir and return it.
fn staged_inputs(fix: &Path, tag: &str) -> PathBuf {
    let tmp = std::env::temp_dir().join(format!("beagle_rs_parity_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(&tmp).unwrap();
    for entry in std::fs::read_dir(fix.join("inputs")).unwrap() {
        let p = entry.unwrap().path();
        std::fs::copy(&p, tmp.join(p.file_name().unwrap())).unwrap();
    }
    tmp
}

/// Run `beagle-rs <scenario args> out=out nthreads=1 seed=99999` in `dir` and assert
/// `out.vcf.gz` is byte-identical to `fixtures/official/golden/<golden>`.
fn assert_beagle_parity(tag: &str, scenario_args: &[&str], golden: &str) {
    let Some(fix) = fixtures_dir() else {
        eprintln!("parity[{tag}]: fixtures absent — skipping");
        return;
    };
    let dir = staged_inputs(&fix, tag);
    let mut args: Vec<String> = scenario_args.iter().map(|s| s.to_string()).collect();
    args.extend(["out=out", "nthreads=1", "seed=99999"].map(String::from));

    let status = Command::new(env!("CARGO_BIN_EXE_beagle-rs"))
        .args(&args)
        .current_dir(&dir)
        .env("SOURCE_DATE_EPOCH", GOLDEN_SOURCE_DATE_EPOCH)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("spawn beagle-rs");
    assert!(
        status.success(),
        "parity[{tag}]: beagle-rs exited with {status}"
    );

    let produced = std::fs::read(dir.join("out.vcf.gz")).expect("out.vcf.gz produced");
    let expected = std::fs::read(fix.join("golden").join(golden)).expect("golden present");
    let ok = produced == expected;
    let _ = std::fs::remove_dir_all(&dir);
    assert!(
        ok,
        "parity[{tag}]: out.vcf.gz ({} bytes) not byte-identical to golden {golden} ({} bytes)",
        produced.len(),
        expected.len()
    );
}

#[test]
#[ignore = "e2e: ~16s debug; run in release via scripts/parity-check.sh"]
fn parity_gt_phasing() {
    assert_beagle_parity("gt", &["gt=test.27Feb25.75f.vcf.gz"], "out.gt.vcf.gz");
}

#[test]
#[ignore = "e2e: run in release via scripts/parity-check.sh"]
fn parity_ref_vcf_imputation() {
    assert_beagle_parity(
        "ref",
        &["ref=ref.27Feb25.75f.vcf.gz", "gt=target.27Feb25.75f.vcf.gz"],
        "out.ref.vcf.gz",
    );
}

#[test]
#[ignore = "e2e: run in release via scripts/parity-check.sh"]
fn parity_ref_bref3_imputation() {
    assert_beagle_parity(
        "bref3imp",
        &["ref=ref.27Feb25.75f.bref3", "gt=target.27Feb25.75f.vcf.gz"],
        "out.bref3.vcf.gz",
    );
}

/// The `bref3` tool must reproduce the Java bref3 encoder byte-for-byte: encoding the reference
/// VCF must equal the committed `inputs/ref.27Feb25.75f.bref3` (itself produced by Java bref3).
#[test]
#[ignore = "e2e: run in release via scripts/parity-check.sh"]
fn parity_bref3_tool() {
    let Some(fix) = fixtures_dir() else {
        eprintln!("parity[bref3tool]: fixtures absent — skipping");
        return;
    };
    let dir = staged_inputs(&fix, "bref3tool");
    let output = Command::new(env!("CARGO_BIN_EXE_bref3"))
        .arg("ref.27Feb25.75f.vcf.gz")
        .current_dir(&dir)
        .stderr(Stdio::null())
        .output()
        .expect("spawn bref3");
    assert!(
        output.status.success(),
        "bref3 exited with {}",
        output.status
    );

    let expected = std::fs::read(fix.join("inputs/ref.27Feb25.75f.bref3")).expect("golden bref3");
    let ok = output.stdout == expected;
    let _ = std::fs::remove_dir_all(&dir);
    assert!(
        ok,
        "parity[bref3tool]: bref3 output ({} bytes) not byte-identical to golden ({} bytes)",
        output.stdout.len(),
        expected.len()
    );
}
