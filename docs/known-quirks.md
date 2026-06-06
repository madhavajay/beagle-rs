# Known Beagle quirks / bugs (preserve-vs-fix log)

These are observed behaviors in the Beagle 5.5 (`27Feb25.75f`) Java source that the
Rust port must **reproduce** (tests are the spec). The preserve-vs-fix decision is
deferred to a single dedicated commit at the end of the port (kestrel-rs pattern).
Default stance: **preserve** (byte-for-byte parity is the goal).

## `ints` package

| ID | Location | Behavior | Decision |
|----|----------|----------|----------|
| ints-1 | `PackedIntArray.get(index)` | Bounds check is `if (index > size)` (uses `>`, not `>=`) and does not check negatives. So `get(size)` does **not** throw and returns whatever the packing formula yields from the backing `int[]`. | preserve |
| ints-2 | `CharArray(byte[] ba)` | Length check is `if ((ba.length % 1) != 0)` — `% 1` is always 0, so the check is dead code (intended `% 2`). Odd-length input then reads `ba[j+1]` out of bounds on the last pair. | preserve |
| ints-3 | `CharArray(int[] ia, int to, int from)` | Parameter names are swapped vs. convention; positionally it is `(start, end)` and computes `new char[end - start]`, looping `start..end`. Functionally normal, but the names mislead. | preserve |
| ints-4 | `IntList.copyOfRange(start, end)` | `new int[end - start]` → `NegativeArraySizeException` when `end < start` (javadoc says "if end > start", a typo). When `end >= size`, copies `size - start` and leaves the tail as 0. | preserve |
| ints-5 | `IntList.get/set(index)` | Only checks `index >= size`, not `index < 0` (negative index → array-index exception instead of the documented `IndexOutOfBoundsException` path). | preserve |
| ints-6 | `PackedIntArray.fromByteArray(ba, from, to, ...)` | Array slot index uses absolute `j` (`ia[j >> indexShift]`) while the bit offset uses relative `offset` (`(offset & valuesPerIntM1)`). Inconsistent when `from != 0`. Only called with `from == 0` in-tree. | preserve |
| ints-7 | `IntIntMap(capacity)` | Code rejects `capacity < 1` (so `0` throws), though javadoc says `capacity < 0`. | preserve |

## `blbutil` package

| ID | Location | Behavior | Decision |
|----|----------|----------|----------|
| blb-1 | `BitArray.getAsInt(index)` | For `index % 64 == 63`, a set bit yields `-1` (the sign bit propagates through the arithmetic `>>`), not `1`. | preserve |
| blb-2 | `BGZIPOutputStream.write(byte[] buf, int off, int len)` | The loop guard `(len - off) >= availSize` and trailing copy of `len` bytes are only correct for `off == 0` with `len` smaller than the block size. Larger or offset writes under-flush and then overflow the fixed `input[]` buffer (`ArrayIndexOutOfBoundsException`). Beagle only ever drives it via small `off==0` writes, so the bug is latent. The Rust port uses correct buffering that flushes at exactly `MAX_INPUT_BYTES`, producing identical block boundaries for all realistic call patterns. | port-correct (latent bug not reproduced; document) |

## `beagleutil` package

| ID | Location | Behavior | Decision |
|----|----------|----------|----------|
| beagleutil-1 | `ChromInterval.isValidPos` | The leading-zero guard `s.charAt(startIndex) == 0` compares to the NUL char (int `0`), not `'0'`, so it never fires — leading zeros in positions are accepted. Preserved. | preserve |

## `vcf` package

| ID | Location | Behavior | Decision |
|----|----------|----------|----------|
| vcf-1 | `AlleleRefGTRec.deepCopy(int[][])` | Builds `copy` but `return ia` (the input), so it is a shallow share, not a deep copy. No observable effect because the record is immutable. The Rust port clones (equivalent for immutable data). | preserve-effect |
| vcf-2 | `LowMafGTRec.alleleCount(majorAllele)` | Subtracts `hapIndices.length` (the allele *count*) once per non-major allele instead of that row's length (`hapIndices[al].length`), so the returned major-allele count is generally wrong. Preserved. | preserve |
| vcf-3 | `PlinkGenMap.closestIndex` | The out-of-range branch tests `insPt == basePos.length` (the chromosome-count / outer-array length) instead of `basePos[chrom].length`, so the right-edge clamp is keyed off the wrong bound. Preserved. | preserve |
