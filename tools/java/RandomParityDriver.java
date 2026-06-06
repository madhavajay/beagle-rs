import blbutil.Utilities;
import java.util.Arrays;
import java.util.Random;

/**
 * Emits a deterministic transcript of java.util.Random outputs and Utilities.shuffle,
 * reproduced by crates/beagle-rs/tests/random_parity.rs. Proves the Rust java.util.Random
 * port matches the JDK bit-for-bit (the determinism foundation for phasing/imputation).
 */
public final class RandomParityDriver {

    public static void main(String[] args) {
        StringBuilder o = new StringBuilder();
        o.append("# random parity (java.util.Random)\n");

        Random r = new Random(42L);

        StringBuilder a = new StringBuilder();
        for (int i = 0; i < 8; ++i) {
            if (i > 0) a.append(',');
            a.append(r.nextInt());
        }
        o.append("nextInt=").append(a).append('\n');

        a = new StringBuilder();
        for (int i = 0; i < 8; ++i) {
            if (i > 0) a.append(',');
            a.append(r.nextInt(100));
        }
        o.append("nextIntBound100=").append(a).append('\n');

        a = new StringBuilder();
        for (int i = 0; i < 8; ++i) {
            if (i > 0) a.append(',');
            a.append(r.nextInt(64)); // power-of-two path
        }
        o.append("nextIntBound64=").append(a).append('\n');

        a = new StringBuilder();
        for (int i = 0; i < 4; ++i) {
            if (i > 0) a.append(',');
            a.append(r.nextLong());
        }
        o.append("nextLong=").append(a).append('\n');

        a = new StringBuilder();
        for (int i = 0; i < 4; ++i) {
            if (i > 0) a.append(',');
            a.append(Double.doubleToLongBits(r.nextDouble())); // exact bits, no float formatting
        }
        o.append("nextDoubleBits=").append(a).append('\n');

        a = new StringBuilder();
        for (int i = 0; i < 8; ++i) {
            if (i > 0) a.append(',');
            a.append(r.nextBoolean() ? 1 : 0);
        }
        o.append("nextBoolean=").append(a).append('\n');

        int[] ia = new int[20];
        for (int i = 0; i < ia.length; ++i) {
            ia[i] = i;
        }
        Utilities.shuffle(ia, new Random(7L));
        o.append("shuffle=").append(Arrays.toString(ia)).append('\n');

        System.out.print(o);
    }
}
