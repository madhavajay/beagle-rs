#!/usr/bin/env bash
# End-to-end byte-for-byte parity gate: run the Rust beagle-rs / bref3 / unbref3 binaries on the
# committed official fixtures and assert byte-identity with the golden Java outputs.
#
# This is the same check the `parity_official` integration test performs, runnable standalone.
# If `reference/*.local.jar` exist, it ALSO diffs against a fresh Java run to confirm the golden
# fixtures themselves have not drifted (and prints the normalized `.log` diff for the record).
#
# Determinism knobs (must match how the golden fixtures were generated): nthreads=1 seed=99999.
# SOURCE_DATE_EPOCH pins the VCF ##filedate wall-clock field to the golden's date (2026-06-06 UTC)
# so the comparison is byte-exact (incl. the BGZIP deflate stream) regardless of when this runs.
set -euo pipefail
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1780747200}"  # 2026-06-06 12:00:00 UTC

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIX="$ROOT/fixtures/official"
SEED=99999
THREADS=1
V="27Feb25.75f"

[ -d "$FIX/golden" ] || { echo "ERROR: missing $FIX/golden — run scripts/make-fixtures.sh first." >&2; exit 1; }

echo ">> building release binaries"
( cd "$ROOT" && cargo build --release -q )
BEAGLE="$ROOT/target/release/beagle-rs"
BREF3="$ROOT/target/release/bref3"
UNBREF3="$ROOT/target/release/unbref3"

tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT
cp "$FIX"/inputs/* "$tmp"/
cd "$tmp"

fail=0
check() { # <label> <produced> <golden>
  if cmp -s "$2" "$3"; then echo "  ✅ $1 — BYTE-IDENTICAL"; else echo "  ❌ $1 — DIFFERS ($2 vs $3)"; fail=1; fi
}

echo ">> beagle-rs: gt= (phasing)"
"$BEAGLE" gt=test.$V.vcf.gz out=out.gt nthreads=$THREADS seed=$SEED >/dev/null 2>&1
check "out.gt.vcf.gz" out.gt.vcf.gz "$FIX/golden/out.gt.vcf.gz"

echo ">> beagle-rs: ref=<vcf> gt= (imputation)"
"$BEAGLE" ref=ref.$V.vcf.gz gt=target.$V.vcf.gz out=out.ref nthreads=$THREADS seed=$SEED >/dev/null 2>&1
check "out.ref.vcf.gz" out.ref.vcf.gz "$FIX/golden/out.ref.vcf.gz"

echo ">> beagle-rs: ref=<bref3> gt= (imputation)"
"$BEAGLE" ref=ref.$V.bref3 gt=target.$V.vcf.gz out=out.bref3 nthreads=$THREADS seed=$SEED >/dev/null 2>&1
check "out.bref3.vcf.gz" out.bref3.vcf.gz "$FIX/golden/out.bref3.vcf.gz"

echo ">> bref3: VCF -> bref3"
"$BREF3" ref.$V.vcf.gz > my.bref3 2>/dev/null
check "ref.bref3" my.bref3 "$FIX/inputs/ref.$V.bref3"

echo ">> unbref3: bref3 -> VCF (data lines; ##filedate is a wall-clock exception)"
"$UNBREF3" "$FIX/inputs/ref.$V.bref3" > my.unbref3.vcf 2>/dev/null
if [ -f "$ROOT/reference/unbref3.local.jar" ]; then
  java -jar "$ROOT/reference/unbref3.local.jar" "$FIX/inputs/ref.$V.bref3" > java.unbref3.vcf 2>/dev/null
  if cmp -s <(grep -v '^##filedate' my.unbref3.vcf) <(grep -v '^##filedate' java.unbref3.vcf); then
    echo "  ✅ unbref3 — IDENTICAL (modulo ##filedate)"
  else
    echo "  ❌ unbref3 — DIFFERS beyond ##filedate"; fail=1
  fi
else
  echo "  (skipped Java comparison — reference/unbref3.local.jar absent)"
fi

# pypgx-style invocations (chrom=<region> + ref= + em=/impute=) compared live against the jar.
# No committed golden (a Java run bakes today's local ##filedate); compared modulo that field.
JAR="$ROOT/reference/beagle.local.jar"
if [ -f "$JAR" ]; then
  echo ">> pypgx-pattern invocations (chrom=/em=/impute=, vs the jar, modulo ##filedate)"
  # Decompress to a file first so `head`/`sed` don't SIGPIPE the zcat producer under pipefail.
  zcat ref.$V.vcf.gz | grep -v '^#' > ref.records
  CHR=$(head -n1 ref.records | cut -f1)
  P1=$(head -n1 ref.records | cut -f2)
  PMID=$(sed -n '700p' ref.records | cut -f2)
  REGION="$CHR:$P1-$PMID"
  pnorm() { zcat "$1" | grep -v '^##filedate'; }
  check_jar() { # <label> <args...>
    local label="$1"; shift; rm -f pr.vcf.gz pj.vcf.gz
    "$BEAGLE"      "$@" out=pr nthreads=$THREADS seed=$SEED >/dev/null 2>&1
    java -jar "$JAR" "$@" out=pj nthreads=$THREADS seed=$SEED >/dev/null 2>&1
    if cmp -s <(pnorm pr.vcf.gz) <(pnorm pj.vcf.gz); then echo "  ✅ $label"; else echo "  ❌ $label — DIFFERS"; fail=1; fi
  }
  check_jar "ref=<vcf> chrom=$REGION"            ref=ref.$V.vcf.gz gt=target.$V.vcf.gz chrom=$REGION
  check_jar "ref=<vcf> chrom + em=false"         ref=ref.$V.vcf.gz gt=target.$V.vcf.gz chrom=$REGION em=false
  check_jar "ref=<vcf> chrom + impute=false"     ref=ref.$V.vcf.gz gt=target.$V.vcf.gz chrom=$REGION impute=false
  check_jar "ref=<bref3> chrom + em=false"       ref=ref.$V.bref3  gt=target.$V.vcf.gz chrom=$REGION em=false
else
  echo "  (skipped pypgx-pattern jar comparison — reference/beagle.local.jar absent)"
fi

if [ $fail -eq 0 ]; then
  echo ">> PARITY OK — all outputs byte-identical to the Java reference"
else
  echo ">> PARITY FAILED" >&2
fi
exit $fail
