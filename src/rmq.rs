#![allow(unused)]
use rayon::iter::{IntoParallelIterator, ParallelIterator};

pub trait Rmq {
    fn name() -> String;
    /// To save time, only run benchmarks up to this n.
    fn max_n() -> usize {
        usize::MAX
    }
    fn build(data: &[u64]) -> Self;
    /// Space usage in bytes.
    fn space(&self) -> usize;
    fn query(&self, data: &[u64], l: usize, r: usize) -> (u64, usize);
}

// -------------------------------------------------------------
// TODO: Implement the Rmq trait for additional data structures.
// -------------------------------------------------------------

/// Block decomposition: sparse table on block minimums, linear scan within blocks.
///
/// For query [l, r]:
/// - Same block: linear scan.
/// - Across blocks: suffix scan of first block, sparse table on interior blocks, prefix scan of last block.
pub struct BlockRmq<const S: usize> {
    block_min_pos: Vec<u8>,
    sparse: SparseTable,
}
impl<const S: usize> Rmq for BlockRmq<S> {
    fn name() -> String {
        format!("BlockRmq<{S}>")
    }
    fn build(data: &[u64]) -> Self {
        assert!(
            S <= 256,
            "Block size S must fit in u8 for position encoding"
        );
        let n = data.len();
        let (block_mins, block_min_pos): (Vec<u64>, Vec<u8>) = (0..(n + S - 1) / S)
            .into_par_iter()
            .map(|b| {
                (b * S..(b * S + S).min(n))
                    .map(|i| (data[i], (i - b * S) as u8))
                    .min()
                    .unwrap()
            })
            .unzip();
        BlockRmq {
            block_min_pos,
            sparse: SparseTable::build(&block_mins),
        }
    }
    fn space(&self) -> usize {
        self.sparse.space()
    }
    fn query(&self, data: &[u64], l: usize, r: usize) -> (u64, usize) {
        let block_l = l / S;
        let block_r = r / S;
        if block_l == block_r {
            return (l..r + 1).map(|i| (data[i], i)).min().unwrap();
        }
        let suffix = (l..(block_l + 1) * S).map(|i| (data[i], i)).min().unwrap();
        let prefix = (block_r * S..r + 1).map(|i| (data[i], i)).min().unwrap();
        let mid = if block_r > block_l + 1 {
            let (val, block_idx) = self.sparse.query(data, block_l + 1, block_r - 1);
            let idx = self.block_min_pos[block_idx];
            (val, block_idx * S + idx as usize)
        } else {
            (u64::MAX, usize::MAX)
        };
        suffix.min(prefix).min(mid)
    }
}

/// Like `BlockRmq`, but precomputes prefix and suffix minima within each block so that
/// cross-block queries need only two O(1) lookups instead of two linear scans.
///
/// `prefix_min[i]` = min(data[block_start ..= i])
/// `suffix_min[i]` = min(data[i ..= block_end])
pub struct BlockRmqPrecomputed<const S: usize> {
    prefix_min: Vec<(u64, usize)>,
    suffix_min: Vec<(u64, usize)>,
    block_min_pos: Vec<usize>,
    sparse: SparseTable,
}
impl<const S: usize> Rmq for BlockRmqPrecomputed<S> {
    fn name() -> String {
        format!("BlockPSRmq<{S}>")
    }
    fn build(data: &[u64]) -> Self {
        let n = data.len();
        let num_blocks = (n + S - 1) / S;
        let mut prefix_min = vec![(0u64, 0); n];
        let mut suffix_min = vec![(0u64, 0); n];

        for b in 0..num_blocks {
            let lo = b * S;
            let hi = (lo + S).min(n);
            // Prefix minima: left to right within block.
            prefix_min[lo] = (data[lo], lo);
            for i in lo + 1..hi {
                prefix_min[i] = prefix_min[i - 1].min((data[i], i));
            }
            // Suffix minima: right to left within block.
            suffix_min[hi - 1] = (data[hi - 1], hi - 1);
            for i in (lo..hi - 1).rev() {
                suffix_min[i] = suffix_min[i + 1].min((data[i], i));
            }
        }

        let (block_mins, block_min_pos): (Vec<u64>, Vec<usize>) = (0..num_blocks)
            .into_par_iter()
            .map(|b| prefix_min[(b * S + S).min(n) - 1])
            .unzip();

        BlockRmqPrecomputed {
            prefix_min,
            suffix_min,
            block_min_pos,
            sparse: SparseTable::build(&block_mins),
        }
    }
    fn space(&self) -> usize {
        std::mem::size_of_val(self.prefix_min.as_slice())
            + std::mem::size_of_val(self.suffix_min.as_slice())
            + self.sparse.space()
    }
    fn query(&self, data: &[u64], l: usize, r: usize) -> (u64, usize) {
        let block_l = l / S;
        let block_r = r / S;
        if block_l == block_r {
            return (l..=r).map(|i| (data[i], i)).min().unwrap();
        }
        let suffix = self.suffix_min[l];
        let prefix = self.prefix_min[r];
        let mid = if block_r > block_l + 1 {
            let (val, block_idx) = self.sparse.query(data, block_l + 1, block_r - 1);
            let idx = self.block_min_pos[block_idx];
            (val, idx)
        } else {
            (u64::MAX, usize::MAX)
        };
        suffix.min(prefix).min(mid)
    }
}

/// Sparse table: O(n log n) space, O(1) query.
///
/// `table[k][i]` = minimum of `data[i .. i + 2^k]`.
/// Query [l, r]: let k = floor(log2(r - l + 1)), return min(table[k][l], table[k][r - 2^k + 1]).
struct SparseTable {
    table: Vec<Vec<(u64, usize)>>,
}
impl Rmq for SparseTable {
    fn name() -> String {
        "SparseTable".to_string()
    }
    fn build(data: &[u64]) -> Self {
        let n = data.len();
        let levels = if n <= 1 {
            1
        } else {
            usize::BITS as usize - n.leading_zeros() as usize
        };

        let mut table: Vec<Vec<(u64, usize)>> = Vec::with_capacity(levels);
        table.push(data.iter().enumerate().map(|(i, x)| (*x, i)).collect());
        for k in 1..levels {
            let half = 1 << (k - 1);
            let prev = &table[k - 1];
            let row = (0..n - (1 << k) + 1)
                .map(|i| prev[i].min(prev[i + half]))
                .collect();
            table.push(row);
        }

        SparseTable { table }
    }
    fn space(&self) -> usize {
        self.table
            .iter()
            .map(|row| std::mem::size_of_val(row.as_slice()))
            .sum::<usize>()
    }
    fn query(&self, _data: &[u64], l: usize, r: usize) -> (u64, usize) {
        let k = usize::BITS as usize - (r - l + 1).leading_zeros() as usize - 1;
        self.table[k][l].min(self.table[k][r - (1 << k) + 1])
    }
}
