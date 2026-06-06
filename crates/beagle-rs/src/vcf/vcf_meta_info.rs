//! Port of `vcf/VcfMetaInfo.java` — a `##key=value` VCF meta-information line.

/// Port of `vcf/VcfMetaInfo.java`.
#[derive(Clone, Debug)]
pub struct VcfMetaInfo {
    line: String,
    key: String,
    value: String,
}

/// `VcfMetaInfo.PREFIX`.
pub const PREFIX: &str = "##";
/// `VcfMetaInfo.DELIMITER`.
pub const DELIMITER: char = '=';

impl VcfMetaInfo {
    /// `new VcfMetaInfo(String line)`.
    pub fn new(line: &str) -> Self {
        let line = line.trim();
        assert!(
            line.starts_with(PREFIX),
            "VCF meta-information line: missing starting \"{PREFIX}\": {line}"
        );
        let index = line.find(DELIMITER).map_or(-1, |i| i as i32);
        assert!(
            index > 0 && index != line.len() as i32 - 1,
            "VCF meta-information line: missing \"{DELIMITER}\""
        );
        let index = index as usize;
        VcfMetaInfo {
            line: line.to_string(),
            key: line[2..index].to_string(),
            value: line[index + 1..].to_string(),
        }
    }

    /// `key()`.
    pub fn key(&self) -> &str {
        &self.key
    }

    /// `value()`.
    pub fn value(&self) -> &str {
        &self.value
    }
}

impl std::fmt::Display for VcfMetaInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.line)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_key_value() {
        let mi = VcfMetaInfo::new("##fileformat=VCFv4.2");
        assert_eq!(mi.key(), "fileformat");
        assert_eq!(mi.value(), "VCFv4.2");
        assert_eq!(mi.to_string(), "##fileformat=VCFv4.2");

        let mi = VcfMetaInfo::new("##INFO=<ID=AC,Number=A,Type=Integer>");
        assert_eq!(mi.key(), "INFO");
        assert_eq!(mi.value(), "<ID=AC,Number=A,Type=Integer>");
    }

    #[test]
    #[should_panic]
    fn rejects_missing_prefix() {
        let _ = VcfMetaInfo::new("fileformat=VCFv4.2");
    }

    #[test]
    #[should_panic]
    fn rejects_missing_delimiter() {
        let _ = VcfMetaInfo::new("##noequals");
    }
}
