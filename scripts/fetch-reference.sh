#!/usr/bin/env bash
# Download the official Beagle 5.5 (27Feb25.75f) binaries + example test data
# used as the differential-testing reference. Everything lands under reference/
# (git-ignored) and is sha256-recorded in reference/SHA256SUMS.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/reference"
BASE="https://faculty.washington.edu/browning/beagle"
mkdir -p "$OUT"

command -v curl >/dev/null 2>&1 || { echo "ERROR: curl not found on PATH." >&2; exit 1; }

# Pinned to Beagle 5.5, version 27Feb25.75f (matches the imported submodule source).
FILES=(
  "beagle.27Feb25.75f.jar"
  "bref3.27Feb25.75f.jar"
  "unbref3.27Feb25.75f.jar"
  "run.beagle.27Feb25.75f.example"
  "beagle_5.5_17Dec24.pdf"
  "test.beagle.vcf.gz"   # bundle the run example unpacks into test/ref/target
)

for f in "${FILES[@]}"; do
  echo ">> fetching $f"
  curl -fsSL -o "$OUT/$f" "$BASE/$f"
done

( cd "$OUT" && sha256sum "${FILES[@]}" > SHA256SUMS )
echo ">> done. artifacts in $OUT:"
ls -la "$OUT"
echo ">> sha256:"; cat "$OUT/SHA256SUMS"
