#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLASSES="$ROOT/reference/classes"; OUT_CLASSES="$ROOT/reference/parity-classes"
FIX_DIR="$ROOT/fixtures/marker"
command -v javac >/dev/null 2>&1 || { echo "ERROR: javac not found." >&2; exit 1; }
[ -d "$CLASSES" ] || { echo "ERROR: $CLASSES missing — run scripts/build-java.sh first." >&2; exit 1; }
mkdir -p "$OUT_CLASSES" "$FIX_DIR"
javac -cp "$CLASSES" -d "$OUT_CLASSES" "$ROOT/tools/java/MarkerParityDriver.java"
java -cp "$CLASSES:$OUT_CLASSES" MarkerParityDriver > "$FIX_DIR/parity.txt"
echo ">> wrote $FIX_DIR/parity.txt ($(wc -l < "$FIX_DIR/parity.txt") lines)"
