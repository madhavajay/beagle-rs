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
