#!/usr/bin/env bash
# Regenerate all cross-language parity fixtures from the Java reference.
# Add a line here per package as its parity driver is created.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
bash "$ROOT/scripts/gen-ints-parity.sh"
bash "$ROOT/scripts/gen-blbutil-parity.sh"
