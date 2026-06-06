//! Port of `bref/BrefWriter.java` — the interface for writing phased, non-missing genotypes
//! to a bref (binary reference format) file. `close()` must be called after the last `write()`
//! to flush buffered data.

use std::rc::Rc;

use crate::vcf::{RefGTRec, Samples};

/// Port of the `bref/BrefWriter.java` interface (`Closeable`).
pub trait BrefWriter {
    /// `samples()`.
    fn samples(&self) -> &Samples;

    /// `write(RefGTRec rec)` — writes the phased genotype data.
    fn write(&mut self, rec: Rc<dyn RefGTRec>);

    /// `close()` — flushes any buffered output and releases resources.
    fn close(&mut self);
}
