# beagle-rs

A Rust port of [Beagle 5.5](https://faculty.washington.edu/browning/beagle/beagle.html)
(`27Feb25.75f`) — genotype phasing and genotype imputation.

**Status:** in progress. The goal is a *complete*, byte-for-byte-exact port of all of
Beagle 5.5 (phasing, imputation, and the `bref3`/`unbref3` tools), verified by a
differential parity suite against the original Java. See [TODO.md](TODO.md) for the full
plan and the explicit Definition of Done.

## Layout

| Path | What |
|---|---|
| `beagle/` | Git submodule: the original Beagle 5.5 Java source (`git@github.com:madhavajay/beagle.git`) — the reference implementation and spec. |
| `scripts/` | Build/fetch helpers for the Java reference. |
| `reference/` | (git-ignored) Built/downloaded Java jars + example data, used for differential testing. |
| `TODO.md` | Phased plan + Definition of Done. |
| `crates/` | The Rust port (added during Phase 5). |

## Licensing

Beagle is GPL v3+ (`beagle/LICENSE`); a `beagle/LICENSE.MIT` notice is retained for the
upstream Broad Institute BGZIP files (not present in this release). This port inherits the
GPL v3 license.

## Building the Java reference

No JDK build tool is required — plain `javac` (Java 8+; tested with JDK 26):

```sh
git submodule update --init --recursive   # if beagle/ is empty
scripts/build-java.sh                      # -> reference/beagle.local.jar
```

To also download the official jars + example test data for differential testing:

```sh
scripts/fetch-reference.sh                 # -> reference/*.jar, test data, SHA256SUMS
```

Run the locally built reference:

```sh
java -jar reference/beagle.local.jar       # prints usage
```

## Building / testing the Rust port

(Added in Phase 5.) Will be the standard:

```sh
cargo build --workspace
cargo test --workspace
```

plus a differential parity gate that runs the Rust binary and the Java reference on the
same inputs and requires byte-for-byte identical output.

## Approach

Test-first and parity-gated: write a high-coverage JUnit suite against the Java source,
capture golden input/output fixtures, then port module-by-module bottom-up
(`ints → blbutil → beagleutil → vcf → bref → main → phase → imp`), mirroring each Java
test in Rust as the spec. Determinism is achievable because Beagle's randomness is plain
`java.util.Random` seeded from `seed=`. Details and gotchas are in [TODO.md](TODO.md).
