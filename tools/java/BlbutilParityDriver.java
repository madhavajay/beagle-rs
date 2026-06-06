import blbutil.BitArray;
import blbutil.StringUtil;
import java.util.Arrays;

/**
 * Emits a deterministic transcript exercising parts of the Beagle {@code blbutil}
 * package whose representation must match byte-for-byte (BitArray word layout) or
 * whose splitting semantics feed VCF parsing (StringUtil). Reproduced by the Rust
 * port in crates/beagle-rs/tests/blbutil_parity.rs.
 */
public final class BlbutilParityDriver {

    public static void main(String[] args) {
        StringBuilder out = new StringBuilder();
        out.append("# blbutil parity transcript (Beagle 5.5 27Feb25.75f)\n");

        // --- BitArray: byte-exact long-word layout + derived ops ---
        BitArray b = new BitArray(200);
        for (int j = 0; j < 200; ++j) {
            if (j % 7 == 0 || j % 11 == 3) {
                b.set(j);
            }
        }
        out.append("bitarray words=").append(Arrays.toString(b.toLongArray())).append('\n');
        out.append("bitarray getAsInt[0,7,63,64,193]=")
           .append(b.getAsInt(0)).append(',')
           .append(b.getAsInt(7)).append(',')
           .append(b.getAsInt(63)).append(',')
           .append(b.getAsInt(64)).append(',')
           .append(b.getAsInt(193)).append('\n');

        BitArray r = b.restrict(5, 130);
        out.append("bitarray restrict5_130 size=").append(r.size())
           .append(" words=").append(Arrays.toString(r.toLongArray())).append('\n');

        out.append("bitarray hash 0_200=").append(b.hash(0, 200))
           .append(" hash 5_130=").append(b.hash(5, 130)).append('\n');

        BitArray d = new BitArray(200);
        d.copyFrom(b, 10, 150);
        out.append("bitarray copyfrom10_150 words=")
           .append(Arrays.toString(d.toLongArray())).append('\n');

        out.append("bitarray longHashCode=")
           .append(BitArray.longHashCode(0x123456789abcdefL)).append('\n');

        // --- StringUtil ---
        out.append("split_tab=").append(join(StringUtil.getFields("a\tb\t\tc\t", '\t'))).append('\n');
        out.append("split_tab_limit2=").append(join(StringUtil.getFields("a\tb\tc\td", '\t', 2))).append('\n');
        out.append("count_tab=").append(StringUtil.countFields("a\tb\tc", '\t'))
           .append(" count_tab_max2=").append(StringUtil.countFields("a\tb\tc\td", '\t', 2)).append('\n');
        out.append("split_ws=").append(join(StringUtil.getFields("  hello   world  foo "))).append('\n');
        out.append("split_ws_limit2=").append(join(StringUtil.getFields("  hello   world  foo ", 2))).append('\n');
        out.append("count_ws=").append(StringUtil.countFields("  hello   world  foo ")).append('\n');

        System.out.print(out);
    }

    private static String join(String[] fields) {
        StringBuilder sb = new StringBuilder();
        for (int j = 0; j < fields.length; ++j) {
            if (j > 0) {
                sb.append('|');
            }
            sb.append(fields[j]);
        }
        return sb.toString();
    }
}
