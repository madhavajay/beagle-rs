# Beagle → Rust Port — TODO

Port **all of** [Beagle 5.5](https://faculty.washington.edu/browning/beagle/beagle.html)
(`27Feb25.75f` — genotype phasing + genotype imputation + the bref3/unbref3 tools) from
Java to Rust.

## Definition of Done (the goal)

This goal is **not complete** until the entire program is ported and verified. No partial
port, no "good enough", no permanent stubs. Concretely, ALL of the following must hold:

- [ ] **Every** Java source file's behavior is reproduced in Rust — all 8 packages
      (`ints blbutil beagleutil vcf bref imp main phase`), all 135 files. Nothing skipped.
- [ ] **Full feature set:** phasing (`gt=`), imputation from VCF ref (`ref=`), imputation
      from bref3 ref, and the standalone `bref3` and `unbref3` tools.
- [ ] **Full CLI surface:** every `key=value` argument, default, and validation/error message
      that Beagle 5.5 accepts behaves identically (see `main/Par.java`).
- [ ] **Byte-for-byte exact output** vs. the Java reference across the entire fixture matrix —
      output `*.vcf.gz`, `*.bref3`, and `*.log` files are *byte-identical* (this includes the
      BGZIP-compressed bytes, not just decompressed content — see [gotchas](#known-gotchas--risks)
      for how we hit byte parity on the deflate stream).
- [ ] **No `todo!`/`unimplemented!`/panicking stubs** anywhere in shipped code paths.
- [ ] The Rust binary is a **drop-in replacement** for `beagle.27Feb25.75f.jar`: same args in →
      same bytes out, same exit codes, same stderr/stdout shape.
- [ ] The full differential **parity suite is 100% green** in CI (not a subset), and the
      mirrored Rust unit tests (ported 1:1 from the Java test suite we write) all pass.
- [ ] `pypgx-rs` integration acceptance: `beagle-rs` replaces the `NotPorted` stub in
      `pypgx-rs/src/external.rs` and pypgx's tests pass against it.

Until every box above is checked, the port is **in progress**. Progress is measured by how
much of this is done, not by whether a milestone "works for a demo".

**Strategy:** *don't port blind.* Build a high-coverage unit-test suite **against the original
Java** first, capture golden input/output fixtures from real runs, get the Java side to ≥90%
line coverage, then port module-by-module to Rust — copying each Java test into Rust as the
spec and gating on byte parity with the Java reference. The phase ordering below is the
**order of work**, not a set of milestones any of which is an acceptable stopping point — the
only acceptable stopping point is the Definition of Done above.

This ordering, the two-crate layout, the parity-gate discipline, and the
"tests are the spec / preserve-bugs-deliberately" mindset are all lifted from our
existing ports — see [Reference ports](#reference-ports--lessons) below.

---

## Upstream source & artifacts

**Project:** Beagle 5.5 — author Brian L. Browning, University of Washington.
**Version imported:** `27Feb25.75f` (source archive `beagle.250227.zip`, released 2025-02-27).
**License:** GNU GPL v3 or later (`LICENSE`); MIT notice retained for Broad Institute
BGZIP files (`LICENSE.MIT`) — note those `net/sf/samtools/` files are **not present** in
this release (BGZIP was reimplemented under GPL in `blbutil/`).
**Landing page (source of truth for all of the below):**
https://faculty.washington.edu/browning/beagle/beagle.html

### Source code
- [x] Java source imported into the `beagle/` submodule (`git@github.com:madhavajay/beagle.git`),
      flattened to top level: 135 `.java` files across 8 packages
      (`beagleutil blbutil bref imp ints main phase vcf`).
- Submodule pinned in this repo's `.gitmodules`.

### Binaries / tools (download for the Java reference + differential testing)
| Artifact | File | URL |
|---|---|---|
| Main program (Java 8) | `beagle.27Feb25.75f.jar` | https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar |
| VCF → bref3 | `bref3.27Feb25.75f.jar` | https://faculty.washington.edu/browning/beagle/bref3.27Feb25.75f.jar |
| bref3 → VCF | `unbref3.27Feb25.75f.jar` | https://faculty.washington.edu/browning/beagle/unbref3.27Feb25.75f.jar |
| Manual (PDF) | `beagle_5.5_17Dec24.pdf` | https://faculty.washington.edu/browning/beagle/beagle_5.5_17Dec24.pdf |
| Source zip | `beagle.250227.zip` | https://faculty.washington.edu/browning/beagle/beagle.250227.zip |
| Run example | `run.beagle.27Feb25.75f.example` | https://faculty.washington.edu/browning/beagle/run.beagle.27Feb25.75f.example |

### Test data (from the official run example)
Downloaded by `run.beagle.27Feb25.75f.example` (bundle at `.../test.beagle.vcf.gz`):
- `test.27Feb25.75f.vcf.gz` — single-panel phasing input
- `ref.27Feb25.75f.vcf.gz` — reference panel (VCF)
- `target.27Feb25.75f.vcf.gz` — target genotypes to impute
- `ref.27Feb25.75f.bref3` — reference panel in bref3 (produced by `bref3.jar`)

Canonical example commands (these become our first golden cases):
```sh
# build bref3 from VCF panel
java -jar bref3.27Feb25.75f.jar ref.27Feb25.75f.vcf.gz > ref.27Feb25.75f.bref3
# phasing only
java -jar beagle.27Feb25.75f.jar gt=test.27Feb25.75f.vcf.gz out=out.gt
# imputation from VCF reference
java -jar beagle.27Feb25.75f.jar ref=ref.27Feb25.75f.vcf.gz gt=target.27Feb25.75f.vcf.gz out=out.ref
# imputation from bref3 reference
java -jar beagle.27Feb25.75f.jar ref=ref.27Feb25.75f.bref3 gt=target.27Feb25.75f.vcf.gz out=out.bref3
```

### Genetic maps & reference panels
- **Genetic maps** (HapMap, GRCh36/37/38): https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/
- **Reference panels** (1000 Genomes phase 3 v5a): https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/

### Secondary data source — `pypgx-rs` (small-scale imputation)
`/home/linux/dev/pypgx-rs` wraps Beagle for pharmacogenomic imputation and is a great
source of *small, real* fixtures.
- Invocation: `pypgx/pypgx/api/utils.py` → `estimate_phase_beagle()` runs
  `java -Xmx2g -jar <jar> gt=<in.vcf> chrom=<region> ref=<panel.vcf.gz> out=<out> impute=<bool> em=<bool>`
  (retries with `em=false` on `IllegalArgumentException: 1.0`).
- Bundled jar: `pypgx/pypgx/api/beagle.22Jul22.46e.jar` — **older** than our imported source
  (`22Jul22` vs `27Feb25.75f`). See [open question](#open-questions) on which version is canonical.
- Reference panel: `~/pypgx-bundle/1kgp/{assembly}/{gene}.vcf.gz` (per-gene 1KGP; bundle downloaded separately).
- Tiny fixtures already in-repo (~700 B each), `VcfFrame[Imported]` zip = `metadata.txt` + `data.vcf`:
  - `pypgx-rs/tests/fixtures/CYP4F2-GRCh37.zip`, `CYP4F2-GRCh38.zip`
  - `pypgx-rs/pypgx/test-data/CYP4F2-GRCh37.zip`, `CYP4F2-GRCh38.zip`
- The Rust side already stubs the call: `pypgx-rs/src/external.rs` returns `PgxError::NotPorted(...)`,
  so a real `beagle-rs` could eventually drop in here.
- Use CYP4F2 as the canonical **minimal** imputation case (region-restricted, 2-variant edge cases exercised).

---

## Reference ports — lessons

Three prior Java/C→Rust ports in `/home/linux/dev/biovault-app/main/repos`. Read their
`TODO.md` before starting; key reusable patterns:

- **`kestrel-rs`** (Java→Rust) — *phased* plan (~14 phases), "tests are the spec":
  every Rust test mirrors a Java test; **preserve documented Java bugs**, decide
  preserve-vs-fix per bug in one dedicated final commit. Two crates (`lib` + caller).
  Coverage gate via `cargo llvm-cov --fail-under-lines`.
- **`samtools-rs`** (C→Rust) — upstream `test.pl` is the *required* parity gate (not a subset).
  Byte-for-byte for binary outputs; semantic parity for headers/stderr/help.
  **Honest-gate discipline:** a missing tool (`bgzip`/`tabix`) silently skipped tests and
  hid 7 real bugs → add a hard preflight check that errors loudly if deps are absent.
- **`bcftools-rs`** (C→Rust) — plugin/module-at-a-time, **one PR per batch, no stacked PRs**,
  each slice green on both Rust + upstream parity gates before merge. Shared engines
  (filter expr) reused across modules. Float rendering must match C `%g` + exact constants.

Distilled conventions to adopt here:
- Two-crate workspace: `crates/beagle-rs` (lib) + `crates/beagle-rs-cli` (thin dispatcher).
- Differential testing against the Java jar is the ultimate gate; unit tests mirror Java tests.
- Fixtures copied in-repo, embedded with `include_bytes!`; deterministic, no network in tests.
- One focused PR per batch; `cargo fmt --check` + `clippy -D warnings` + `cargo test` + coverage all green.
- Keep this TODO in sync with CI after every landed batch.

---

## Phase 0 — Repo & harness setup  ✅ DONE (2026-06-06)
Toolchain present (2026-06-06): **Java 26** (`javac` 26 — compiles the Java 8 source fine), **Rust 1.95**,
`curl` (no `wget`), **no** `mvn`/`gradle`/`ant`, **no** `bgzip`/`tabix` (Beagle ships its own BGZIP — not needed).
- [x] `README.md` — project overview, layout, build/test instructions.
- [x] `scripts/fetch-reference.sh` — `curl` down `beagle/bref3/unbref3` jars + example + test data into
      `reference/` (git-ignored), with `SHA256SUMS`.
- [x] `scripts/build-java.sh` — compiles `beagle/**.java` to `reference/beagle.local.jar` with plain `javac`.
- [x] Preflight checks: both scripts hard-error if `javac`/`curl` missing (samtools-rs lesson).

## Phase 1 — Build & run the Java reference  ✅ DONE (2026-06-06)
- [x] Compiled `beagle/` (135 files → 148 classes, 0 errors) to `reference/beagle.local.jar`.
- [x] **Verified the local build is byte-for-byte identical to the official jar** across all three modes
      (`gt=` phasing, `ref=` VCF imputation, bref3 imputation) — confirms imported source = the exact 5.5
      release. Local `bref3` tool also produces a **byte-identical `.bref3`** (62481 B) vs the official tool.
- [x] **Determinism confirmed:** with `nthreads=1 seed=99999`, output is identical across repeated runs →
      byte-for-byte parity testing is achievable. (Default `seed=99999`; default `nthreads`=cores.)
- [x] Characterized `unbref3` round-trip: bref3 deliberately drops VCF header metadata/INFO and regenerates
      a canonical minimal header (`VCFv4.2` + fresh `filedate`); haplotype/genotype data is lossless. The
      Rust `unbref3` must reproduce this canonicalization, not the original header.
- Note: default multi-thread runs not yet checked for determinism — port targets `nthreads=1` first.

## Phase 2 — Test data acquisition  (official ✅ / pypgx pending)
- [x] Vendored the official tiny example data (`test/ref/target` + `ref.bref3`) under `fixtures/official/inputs/`.
- [x] `scripts/make-fixtures.sh` runs the Java jar with fixed `seed`/`nthreads` and stores golden
      `out.{gt,ref,bref3}.vcf.gz` under `fixtures/official/golden/` + `MANIFEST.md` (commands + sha256).
- [x] Characterized the bref3 round-trip (see Phase 1) — `unbref3` canonicalizes the VCF header.
- [ ] Extract `CYP4F2` VCFs from the pypgx zips into `fixtures/pypgx/`; obtain/trim a matching 1KGP
      micro-panel so the run is self-contained and fast. **Blocked on the pypgx-bundle 1KGP panel**
      (downloaded separately; not present locally) — needed to produce golden imputation output.

## Phase 3 — Java unit tests → ~90% line coverage  ← gate before any Rust
The Java source ships **no tests** (confirmed: 0 `*Test*` files, no JUnit). Build them.
- [ ] Add JUnit 5 + JaCoCo on a branch of the `madhavajay/beagle` fork (the submodule repo) so tests
      version alongside the code they cover ([decided](#decisions-resolved-2026-06-06)).
- [ ] Write unit tests bottom-up, package by package (see port order below). Each test documents the
      *observed* behavior — including any bugs — and is the spec the Rust port must reproduce.
- [ ] Track coverage with JaCoCo; **target ≥ 90% line coverage** overall before porting begins.
- [ ] Maintain a **Known-bugs / quirks** log (preserve-vs-fix decided later, per kestrel-rs).

## Phase 4 — Golden / differential harness
- [ ] CLI parity runner: invoke Java jar and (later) Rust binary on identical inputs with fixed
      `seed=` and `nthreads=1`, then compare outputs.
- [ ] **Gate = byte-for-byte identical files** (`*.vcf.gz`, `*.bref3`). Hitting this on `.vcf.gz`
      requires replicating Java's DEFLATE + Beagle's BGZIP block layout exactly — see
      [gotchas](#known-gotchas--risks).
- [ ] Keep a decompressed-content diff as a **diagnostic** to localize mismatches — it does NOT lower
      the bar; the compressed bytes must match too.
- [ ] `*.log`: byte-identical except for intrinsically wall-clock fields (run date, elapsed time). Those
      are the *only* permitted normalization, and each normalized field must be explicitly enumerated
      and justified — nothing else gets a pass.

## Phase 5 — Rust workspace scaffolding  ✅ DONE (2026-06-06)
- [x] `cargo` workspace: `crates/beagle-rs` (lib) + `crates/beagle-rs-cli` (bin, placeholder CLI).
- [x] CI (`.github/workflows/ci.yml`): **rust** gate (`fmt --check`, `clippy -D warnings`, `test`) +
      **java-parity** gate (rebuild Java reference from the submodule, regenerate the parity fixtures,
      fail on drift). Both required.
- [ ] Deps added per-package as needed (e.g. zlib-backed `flate2` for `blbutil` BGZIP); `cargo llvm-cov`
      coverage gate to be added once there is more surface to measure.

### Per-package parity harness pattern (established with `ints`)
For each pure-logic package: a Java driver under `tools/java/` emits a deterministic transcript from the
original classes (`scripts/gen-*-parity.sh` → `fixtures/<pkg>/parity.txt`); a Rust integration test
reproduces it via `include_str!` and asserts byte-for-byte equality. The `ints` transcript includes the
**byte-exact `PackedIntArray` backing words** and the **internal `IntIntMap` key/value order**, so it
verifies representation, not just observable outputs.

## Phase 6 — Port module-by-module (bottom-up)
Port order follows the dependency graph; copy the matching Java tests into Rust first, then implement
until they pass, then add the parity check. File counts in parens.
- [x] `ints` (9) ✅ — `IntArray` trait + `IntList`, `SynchedIntList`, `PackedIntArray`, `UnsignedByteArray`,
      `CharArray`, `WrappedIntArray`, `IndexArray`, `IntIntMap` + factories/statics. 32 unit tests +
      cross-language parity (byte-exact packing + map ordering). 7 quirks preserved (see `docs/known-quirks.md`).
- [ ] `blbutil` (22) — base utils: `BitArray`, `FloatArray/DoubleArray`, `StringUtil`, `Validate`,
      file IO, **BGZIP** (`BGZipIt`, `BGZIPOutputStream`), `MultiThreadUtils`.
- [ ] `beagleutil` (8) — `ChromIds`, `SampleIds`, `ThreadSafeIndexer`, PBWT updaters (`PbwtUpdater`,
      `PbwtDivUpdater`), intervals.
- [ ] `vcf` (38) — VCF/genotype model + IO: `GT`, `Marker(s)`, `VcfRec`, `VcfRecGTParser`, `RefGT*`,
      sliding windows, `PlinkGenMap`, `MarkerMap`. Largest IO surface; many record encodings.
- [ ] `bref` (10) — bref3 binary format read/write (`Bref3*`, `SeqCoder3`, `AsIs/CompressBref3Writer`,
      `UnBref3`). **Byte-exact** target.
- [ ] `main` (5) — `Par` (CLI arg parsing — define the full arg surface), `Main`, `Pedigree`,
      `RunStats`, `WindowWriter`. Match arg names/defaults exactly.
- [ ] `phase` (30) — phasing engine: PBWT/IBS (`Pbwt*Phaser`, `Ibs2*`), HMM (`PhaseBaum*`,
      `HmmUpdater`), Li-Stephens (`PhaseLS`), `SamplePhase`, `FixedPhaseData`. Core algorithm; FP-sensitive.
- [ ] `imp` (12) — imputation engine: `ImpData`, `ImpIbs`, `ImpLS`/`ImpLSBaum`, `HaplotypeCoder`,
      `ImputedRecBuilder`, `ImputedVcfWriter`, `StateProbs`. Core algorithm; FP-sensitive.

## Phase 7 — End-to-end parity & integration
- [ ] Full-pipeline parity on every fixture (official + pypgx CYP4F2): phasing, VCF-ref imputation,
      bref3-ref imputation, bref3 round-trip.
- [ ] Resolve the **preserve-vs-fix** bug decisions in one dedicated commit (kestrel-rs pattern).
- [ ] (Optional) Wire `beagle-rs` into `pypgx-rs/src/external.rs` to replace the `NotPorted` stub and
      re-run pypgx's tests as an external acceptance gate.

---

## Known gotchas & risks
- **Randomness/seed** — Beagle uses `java.util.Random` *only* (confirmed: no `ThreadLocalRandom`/
  `SplittableRandom`/`Math.random`), always re-seeded deterministically from `par.seed()` with fixed
  offsets (`seed + sample`, `seed + step`, `new Random(hap)`, ...). Implications: (1) replicate
  `java.util.Random`'s 48-bit LCG byte-for-byte in Rust (fully specified — `nextInt`/`nextLong`/
  `nextDouble`/bounded rejection sampling are all documented); (2) RNG draws are
  **thread-count-independent**. Pin `seed=` on every reference run.
- **Threading nondeterminism** — `nthreads=` won't change *which* random numbers are drawn (see above)
  but can still reorder floating-point accumulation. Use `nthreads=1` for golden runs; only parallelize
  the Rust port once single-thread parity holds.
- **Floating-point accumulation order** — HMM / Li-Stephens forward-backward sums must match Java's
  order and types (float vs double). Mirror Java arithmetic before "optimizing".
- **bref3 binary format** — custom bit-packed/compressed layout; endianness and field order must be
  byte-exact. Round-trip test early.
- **BGZIP / byte-identical compression** — to hit byte-for-byte parity on `*.vcf.gz` we must reproduce
  `blbutil/BGZIPOutputStream` exactly: identical DEFLATE output (Java `Deflater` is zlib → use a
  zlib-backed deflate in Rust, e.g. `flate2` with the `zlib`/`zlib-ng` backend, **not** pure-Rust
  miniz_oxide, which emits different bytes), same compression level/strategy, same uncompressed block
  size, same BGZF EOF marker, and `MTIME=0` in every block header. Verify the block layout against the
  Java writer byte-by-byte early — this is the single biggest threat to the byte-for-byte goal.
- **VCF record encodings** — many specialized `RefGTRec` variants (allele/bitset/low-MAF/two-allele);
  cover each encoding's parse/round-trip.
- **Java numeric semantics** — signed `int`/overflow, `>>>` vs `>>`, integer division, `String.hashCode`,
  `Float.toString`/`%g`-style formatting (samtools-rs float-formatting lesson).
- **CLI arg surface** — `key=value` style (`gt=`, `ref=`, `map=`, `out=`, `chrom=`, `impute=`, `em=`,
  `nthreads=`, `seed=`, `window=`, ...). Match names, defaults, and validation messages.

## Decisions (resolved 2026-06-06)
- **Canonical version:** `27Feb25.75f` (Beagle 5.5) — the imported source. Regenerate fresh golden
  outputs for the pypgx CYP4F2 data through 5.5 rather than trusting pypgx's older `22Jul22.46e` jar.
- **Java unit tests:** on a branch of the `madhavajay/beagle` fork (the submodule repo), versioned
  alongside the code they cover.
- **Scope:** **full Beagle** — phasing + imputation + bref3/unbref3 tools.
- **RNG:** confirmed `java.util.Random` only, deterministically seeded from `seed=`; replicable in Rust
  and thread-count-independent (see [gotchas](#known-gotchas--risks)).

## Open questions
- [ ] How small can the 1KGP micro-panel be trimmed while still exercising imputation realistically?
      (Empirical — settle during Phase 2.)
