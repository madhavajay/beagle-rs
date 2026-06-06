//! Port of `beagleutil/ChromIds.java` — the process-global registry of chromosome
//! identifiers. The Java `ConcurrentMap`/`volatile List` fast-path caches are perf
//! optimizations over the indexer and are omitted (the indexer is the source of truth).

use super::ThreadSafeIndexer;
use std::sync::OnceLock;

static CHROM_IDS: OnceLock<ChromIds> = OnceLock::new();

/// Port of `beagleutil/ChromIds.java` (singleton).
pub struct ChromIds {
    indexer: ThreadSafeIndexer<String>,
}

impl ChromIds {
    fn create() -> Self {
        ChromIds {
            indexer: ThreadSafeIndexer::with_capacity(4),
        }
    }

    /// `ChromIds.instance()`.
    pub fn instance() -> &'static ChromIds {
        CHROM_IDS.get_or_init(ChromIds::create)
    }

    /// `getIndex(String id)`.
    pub fn get_index(&self, id: &str) -> i32 {
        assert!(!id.is_empty(), "id.isEmpty()");
        self.indexer.get_index(id.to_string())
    }

    /// `getIndices(String[] ids)`.
    pub fn get_indices(&self, ids: &[String]) -> Vec<i32> {
        for id in ids {
            assert!(!id.is_empty(), "id.isEmpty()");
        }
        self.indexer.get_indices(ids)
    }

    /// `getIndexIfIndexed(String id)`.
    pub fn get_index_if_indexed(&self, id: &str) -> i32 {
        assert!(!id.is_empty(), "id.isEmpty()");
        self.indexer.get_index_if_indexed(&id.to_string())
    }

    /// `size()`.
    pub fn size(&self) -> i32 {
        self.indexer.size()
    }

    /// `id(int index)`.
    pub fn id(&self, index: i32) -> String {
        self.indexer.item(index)
    }

    /// `ids()`.
    pub fn ids(&self) -> Vec<String> {
        self.indexer.items()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_chrom_ids() {
        // Uses unique ids since the registry is a process-global singleton.
        let c = ChromIds::instance();
        let i = c.get_index("chrTEST_A");
        assert_eq!(c.id(i), "chrTEST_A");
        assert_eq!(c.get_index("chrTEST_A"), i); // stable
        assert_eq!(c.get_index_if_indexed("chrTEST_A"), i);
        assert_eq!(c.get_index_if_indexed("chrNEVER_SEEN_X"), -1);
    }

    #[test]
    #[should_panic]
    fn rejects_empty_id() {
        let _ = ChromIds::instance().get_index("");
    }
}
