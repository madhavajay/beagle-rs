//! Port of `beagleutil/SampleIds.java` — the process-global registry of sample
//! identifiers.

use super::ThreadSafeIndexer;
use std::sync::OnceLock;

static SAMPLE_IDS: OnceLock<SampleIds> = OnceLock::new();

/// Port of `beagleutil/SampleIds.java` (singleton).
pub struct SampleIds {
    indexer: ThreadSafeIndexer<String>,
}

impl SampleIds {
    fn create() -> Self {
        SampleIds {
            indexer: ThreadSafeIndexer::with_capacity(5000),
        }
    }

    /// `SampleIds.instance()`.
    pub fn instance() -> &'static SampleIds {
        SAMPLE_IDS.get_or_init(SampleIds::create)
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

    /// `idIndexToIndex(Samples samples)` — maps sample-id-index → sample-index, `-1`
    /// for ids absent from `ids`. Takes the sample ids directly (what `Samples.ids()`
    /// returns) to avoid depending on the not-yet-ported `vcf::Samples`.
    pub fn id_index_to_index(&self, ids: &[String]) -> Vec<i32> {
        let id_index: Vec<i32> = ids.iter().map(|id| self.get_index(id)).collect();
        let mut id_index_to_index = vec![-1; self.size() as usize];
        for (j, &ii) in id_index.iter().enumerate() {
            if id_index_to_index[ii as usize] != -1 {
                panic!("Dupicate sample: {}", ids[j]);
            }
            id_index_to_index[ii as usize] = j as i32;
        }
        id_index_to_index
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

    /// `ids(int[] indices)`.
    pub fn ids_at(&self, indices: &[i32]) -> Vec<String> {
        self.indexer.items_at(indices)
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
    fn round_trips_sample_ids() {
        let s = SampleIds::instance();
        let i = s.get_index("sampTEST_1");
        assert_eq!(s.id(i), "sampTEST_1");
        assert_eq!(s.get_index("sampTEST_1"), i);
        let idxs = s.get_indices(&["sampTEST_2".into(), "sampTEST_3".into()]);
        assert_eq!(
            s.ids_at(&idxs),
            vec!["sampTEST_2".to_string(), "sampTEST_3".to_string()]
        );
    }
}
