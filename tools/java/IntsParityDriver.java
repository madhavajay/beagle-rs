import ints.CharArray;
import ints.IntArray;
import ints.IntIntMap;
import ints.IntList;
import ints.PackedIntArray;
import ints.UnsignedByteArray;
import ints.WrappedIntArray;
import java.lang.reflect.Field;
import java.util.Arrays;

/**
 * Emits a deterministic transcript exercising the Beagle {@code ints} package.
 * The Rust port reproduces this transcript byte-for-byte in
 * crates/beagle-rs/tests/ints_parity.rs. Run against reference/classes.
 */
public final class IntsParityDriver {

    private static final int[] VALUE_SIZES =
            {1, 2, 3, 4, 5, 8, 15, 16, 17, 100, 256, 1000, 65536, 70000};

    public static void main(String[] args) throws Exception {
        StringBuilder out = new StringBuilder();
        out.append("# ints parity transcript (Beagle 5.5 27Feb25.75f)\n");

        // --- IntArray.packedCreate / create: get() outputs round-trip the inputs ---
        for (int vs : VALUE_SIZES) {
            int[] values = values(vs);
            IntArray packed = IntArray.packedCreate(values, vs);
            out.append("packedCreate vs=").append(vs)
               .append(" get=").append(getCsv(packed)).append('\n');
            IntArray created = IntArray.create(values, vs);
            out.append("create vs=").append(vs)
               .append(" get=").append(getCsv(created)).append('\n');
        }

        // --- PackedIntArray backing words (byte-exact bit packing) ---
        Field iaField = PackedIntArray.class.getDeclaredField("ia");
        iaField.setAccessible(true);
        for (int vs : VALUE_SIZES) {
            if (vs <= 16) {
                PackedIntArray p = new PackedIntArray(values(vs), vs);
                int[] backing = (int[]) iaField.get(p);
                out.append("packedwords vs=").append(vs)
                   .append(" ia=").append(Arrays.toString(backing)).append('\n');
            }
        }

        // --- IntList ---
        IntList il = new IntList();
        for (int j = 0; j < 30; ++j) {
            il.add(j);
        }
        out.append("intlist toString=").append(il.toString()).append('\n');
        out.append("intlist copyOf35=").append(Arrays.toString(il.copyOf(35))).append('\n');
        out.append("intlist copyOfRange5_40=")
           .append(Arrays.toString(il.copyOfRange(5, 40))).append('\n');
        il.sort();
        out.append("intlist sorted=").append(il.toString()).append('\n');
        out.append("intlist bsearch13=").append(il.binarySearch(13))
           .append(" bsearch1000=").append(il.binarySearch(1000)).append('\n');

        // --- IntIntMap: keys()/values() order is internal -> strong parity check ---
        IntIntMap m = new IntIntMap(4);
        for (int i = 0; i < 300; ++i) {
            int key = (i * 37) % 50;
            if (i % 5 == 0) {
                m.remove(key);
            } else {
                m.put(key, i);
            }
        }
        out.append("intintmap size=").append(m.size())
           .append(" keys=").append(Arrays.toString(m.keys()))
           .append(" values=").append(Arrays.toString(m.values())).append('\n');

        // Demonstrate the UnsignedByteArray/CharArray/WrappedIntArray direct paths too.
        out.append("ubyte=").append(getCsv(new UnsignedByteArray(new int[] {0, 127, 255})))
           .append('\n');
        out.append("char=").append(getCsv(new CharArray(new int[] {0, 4660, 65535})))
           .append('\n');
        out.append("wrapped=").append(getCsv(new WrappedIntArray(new int[] {7, 8, 9})))
           .append('\n');

        System.out.print(out);
    }

    private static int[] values(int valueSize) {
        int[] values = new int[50];
        for (int j = 0; j < values.length; ++j) {
            values[j] = j % valueSize;
        }
        return values;
    }

    private static String getCsv(IntArray a) {
        StringBuilder sb = new StringBuilder();
        for (int j = 0, n = a.size(); j < n; ++j) {
            if (j > 0) {
                sb.append(',');
            }
            sb.append(a.get(j));
        }
        return sb.toString();
    }
}
