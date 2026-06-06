//! Port of `blbutil/SampleFileIt.java` and `blbutil/VcfFileIt.java` (interfaces) — file
//! iterators that also expose the sample list / VCF header.

use super::FileIt;
use crate::vcf::{Samples, VcfHeader};

/// Port of `blbutil/SampleFileIt.java` — a `FileIt` whose records share a sample list.
pub trait SampleFileIt: FileIt {
    /// `samples()`.
    fn samples(&self) -> &Samples;
}

/// Port of `blbutil/VcfFileIt.java` — a `SampleFileIt` that also exposes the VCF header.
pub trait VcfFileIt: SampleFileIt {
    /// `vcfHeader()`.
    fn vcf_header(&self) -> &VcfHeader;
}
