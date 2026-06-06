#!/usr/bin/env bash
# Build the Beagle Java reference from the `beagle/` submodule source using plain
# javac (no build tool is installed). Produces reference/beagle.local.jar, a
# runnable jar equivalent to the published beagle.27Feb25.75f.jar.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$ROOT/beagle"
OUT="$ROOT/reference"
CLASSES="$OUT/classes"
JAR="$OUT/beagle.local.jar"

# --- preflight: fail loudly, never silently skip (samtools-rs lesson) ---
command -v javac >/dev/null 2>&1 || { echo "ERROR: javac not found on PATH (need a JDK)." >&2; exit 1; }
[ -d "$SRC" ] && [ -n "$(find "$SRC" -name '*.java' -print -quit)" ] || {
  echo "ERROR: no Java sources under $SRC — is the 'beagle' submodule checked out?" >&2
  echo "       run: git submodule update --init --recursive" >&2
  exit 1
}

echo ">> javac: $(javac -version 2>&1)"
rm -rf "$CLASSES" && mkdir -p "$CLASSES"

SOURCES="$(mktemp)"
trap 'rm -f "$SOURCES"' EXIT
find "$SRC" -name '*.java' > "$SOURCES"
echo ">> compiling $(wc -l < "$SOURCES" | tr -d ' ') source files"
javac -d "$CLASSES" @"$SOURCES"

echo ">> packaging $JAR (main-class main.Main)"
( cd "$CLASSES" && jar --create --file "$JAR" --main-class main.Main . )

echo ">> built: $JAR"
java -jar "$JAR" 2>&1 | head -1 || true
