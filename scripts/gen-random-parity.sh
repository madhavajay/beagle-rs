#!/usr/bin/env bash
# Compile + run the Java Random parity driver against the reference classes.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLASSES="$ROOT/reference/classes"
DRIVER="$ROOT/tools/java/RandomParityDriver.java"
FIX_DIR="$ROOT/fixtures/random"
OUT_CLASSES="$ROOT/reference/parity-classes"
command -v javac >/dev/null 2>&1 || { echo "ERROR: javac not found." >&2; exit 1; }
[ -d "$CLASSES" ] || { echo "ERROR: $CLASSES missing — run scripts/build-java.sh first." >&2; exit 1; }
mkdir -p "$OUT_CLASSES" "$FIX_DIR"
javac -cp "$CLASSES" -d "$OUT_CLASSES" "$DRIVER"
java -cp "$CLASSES:$OUT_CLASSES" RandomParityDriver > "$FIX_DIR/parity.txt"
echo ">> wrote $FIX_DIR/parity.txt ($(wc -l < "$FIX_DIR/parity.txt") lines)"
