#!/usr/bin/env bash
# Compile + run the Java BGZIP parity driver to produce a reference .bgz the Rust port
# must reproduce byte-for-byte. Requires scripts/build-java.sh first.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLASSES="$ROOT/reference/classes"
DRIVER="$ROOT/tools/java/BgzipParityDriver.java"
FIX_DIR="$ROOT/fixtures/bgzip"
OUT_CLASSES="$ROOT/reference/parity-classes"

command -v javac >/dev/null 2>&1 || { echo "ERROR: javac not found." >&2; exit 1; }
[ -d "$CLASSES" ] || { echo "ERROR: $CLASSES missing — run scripts/build-java.sh first." >&2; exit 1; }

mkdir -p "$OUT_CLASSES" "$FIX_DIR"
javac -cp "$CLASSES" -d "$OUT_CLASSES" "$DRIVER"
java -cp "$CLASSES:$OUT_CLASSES" BgzipParityDriver "$FIX_DIR/sample.bgz"

echo ">> wrote $FIX_DIR/sample.bgz"
echo "   $(wc -c < "$FIX_DIR/sample.bgz") bytes, sha256 $(sha256sum "$FIX_DIR/sample.bgz" | cut -d' ' -f1)"
