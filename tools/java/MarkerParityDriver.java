import vcf.Marker;
import vcf.MarkerParser;

/**
 * Parses a set of VCF records through Beagle's {@code Marker} and emits a transcript of
 * every accessor, for two MarkerParser configurations (store-all and store-none).
 * Reproduced by crates/beagle-rs/tests/marker_parity.rs. Tabs in alleles/toString are
 * rendered as '/' to keep one record per line.
 */
public final class MarkerParityDriver {

    private static final String[] RECORDS = {
        "chr1\t100\trs1\tA\tC\t30\tPASS\t.\tGT\t0|1",
        "chr1\t200\t.\tA\tG,T\t.\t.\tAC=2;AN=4\tGT\t1|2",
        "chr1\t300\tidDEL\tAC\tA\t.\tq10\tEND=305;SVTYPE=DEL\tGT\t0|0",
        "chr1\t400\t.\tA\t.\t.\t.\t.\tGT\t0|0",
        "chr2\t500\trs5\tACGT\tA,ATTT\t99\tPASS\tDP=10\tGT\t0|1",
        "chr2\t600\t.\tG\tA,C,T\t.\tPASS\t.\tGT\t0|3",
    };

    public static void main(String[] args) {
        StringBuilder o = new StringBuilder();
        o.append("# marker parity (Beagle 5.5 27Feb25.75f)\n");
        emit(o, "full", new MarkerParser(true, true, true, true));
        emit(o, "none", new MarkerParser(false, false, false, false));
        System.out.print(o);
    }

    private static void emit(StringBuilder o, String mode, MarkerParser parser) {
        for (int i = 0; i < RECORDS.length; ++i) {
            Marker m = Marker.instance(RECORDS[i], parser);
            o.append("mode=").append(mode).append(" rec=").append(i)
             .append(" chrom=").append(m.chrom())
             .append(" pos=").append(m.pos())
             .append(" id=").append(m.id())
             .append(" alleles=").append(clean(m.alleles()))
             .append(" nAlleles=").append(m.nAlleles())
             .append(" nRefBases=").append(m.nRefBases())
             .append(" qual=").append(m.qual())
             .append(" filter=").append(m.filter())
             .append(" info=").append(m.info())
             .append(" end=").append(m.endValue())
             .append(" hasId=").append(m.hasIdData())
             .append(" hasEnd=").append(m.hasEndValue())
             .append(" bits=").append(m.bitsPerAllele())
             .append(" str=").append(clean(m.toString()))
             .append('\n');
        }
    }

    private static String clean(String s) {
        return s.replace('\t', '/');
    }
}
