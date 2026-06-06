import blbutil.BGZIPOutputStream;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Writes BGZF(deterministic input) using Beagle's {@code BGZIPOutputStream} to the file
 * named by args[0]. The Rust port compresses the same input and compares byte-for-byte
 * (crates/beagle-rs/tests/bgzip_parity.rs). Input is written in small chunks, matching
 * how Beagle's buffered writers actually call the stream.
 */
public final class BgzipParityDriver {

    public static void main(String[] args) throws IOException {
        int n = 200_000;
        byte[] input = new byte[n];
        for (int i = 0; i < n; ++i) {
            input[i] = (byte) ((i * 31 + 7) % 256);
        }
        try (FileOutputStream fos = new FileOutputStream(args[0]);
                BufferedOutputStream bos = new BufferedOutputStream(fos);
                BGZIPOutputStream os = new BGZIPOutputStream(bos, true)) {
            // Write byte-by-byte: the well-defined path that fills each block to exactly
            // MAX_INPUT_BYTES and flushes, giving deterministic 65505-byte block
            // boundaries on both sides. (Beagle's write(buf,off,len) is only safe for
            // small off=0 writes; this isolates the deflate-body parity question.)
            for (int i = 0; i < n; ++i) {
                os.write(input[i]);
            }
        }
    }
}
