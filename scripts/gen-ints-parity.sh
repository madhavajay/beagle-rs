#!/usr/bin/env bash
# Compile + run the Java ints parity driver against the reference classes and write
# the transcript fixture that the Rust port must reproduce byte-for-byte.
# Requires scripts/build-java.sh to have been run (provides reference/classes).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLASSES="$ROOT/reference/classes"
DRIVER="$ROOT/tools/java/IntsParityDriver.java"
FIX_DIR="$ROOT/fixtures/ints"
OUT_CLASSES="$ROOT/reference/parity-classes"

command -v javac >/dev/null 2>&1 || { echo "ERROR: javac not found." >&2; exit 1; }
[ -d "$CLASSES" ] || { echo "ERROR: $CLASSES missing — run scripts/build-java.sh first." >&2; exit 1; }

mkdir -p "$OUT_CLASSES" "$FIX_DIR"
javac -cp "$CLASSES" -d "$OUT_CLASSES" "$DRIVER"
java -cp "$CLASSES:$OUT_CLASSES" IntsParityDriver > "$FIX_DIR/parity.txt"

echo ">> wrote $FIX_DIR/parity.txt"
echo "   $(wc -l < "$FIX_DIR/parity.txt") lines, sha256 $(sha256sum "$FIX_DIR/parity.txt" | cut -d' ' -f1)"
