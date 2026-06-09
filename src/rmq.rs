#![allow(unused)]
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::{lcp::CompactLcp, SaElem, SA};

pub trait RmqElem: Copy + Ord + std::fmt::Debug + Send + Sync {
    const MAX: Self;
}

impl RmqElem for u16 {
    const MAX: Self = u16::MAX;
}

impl RmqElem for u32 {
    const MAX: Self = u32::MAX;
}

impl RmqElem for u64 {
    const MAX: Self = u64::MAX;
}

impl RmqElem for usize {
    const MAX: Self = usize::MAX;
}

pub trait Data<T>: Send + Sync + Copy {
    fn lenn(self) -> usize;
    fn index(self, idx: usize) -> T;
}
impl<T: Copy + Sync + Send> Data<T> for &[T] {
    fn lenn(self) -> usize {
        self.len()
    }

    fn index(self, idx: usize) -> T {
        self[idx]
    }
}
impl Data<u32> for (&SA, &CompactLcp) {
    fn lenn(self) -> usize {
        self.0.len()
    }

    fn index(self, idx: usize) -> u32 {
        self.1.get(&self.0, idx)
    }
}

pub trait Rmq<T: RmqElem> {
    fn name() -> String;
    /// To save time, only run benchmarks up to this n.
    fn max_n() -> usize {
        usize::MAX
    }
    fn build(data: impl Data<T>) -> Self;
    /// Space usage in bytes.
    fn space(&self) -> usize;
    fn query(&self, data: impl Data<T>, l: usize, r: usize) -> (T, usize);
}

// -------------------------------------------------------------
// TODO: Implement the Rmq trait for additional data structures.
// -------------------------------------------------------------

/// Block decomposition: sparse table on block minimums, linear scan within blocks.
///
/// For query [l, r]:
/// - Same block: linear scan.
/// - Across blocks: suffix scan of first block, sparse table on interior blocks, prefix scan of last block.
pub struct BlockRmq<T: RmqElem, const S: usize> {
    block_min_pos: Vec<u8>,
    // sparse: SparseTable<T>,
    segtree: SegmentTree<T>,
}
impl<T: RmqElem, const S: usize> Rmq<T> for BlockRmq<T, S> {
    fn name() -> String {
        format!("BlockRmq<{S}>")
    }
    fn build(data: impl Data<T>) -> Self {
        assert!(
            S <= 256,
            "Block size S must fit in u8 for position encoding"
        );
        let n = data.lenn();
        let (block_mins, block_min_pos): (Vec<T>, Vec<u8>) = (0..(n + S - 1) / S)
            .into_par_iter()
            .map(|b| {
                (b * S..(b * S + S).min(n))
                    .map(|i| (data.index(i), (i - b * S) as u8))
                    .min()
                    .unwrap()
            })
            .unzip();
        BlockRmq {
            block_min_pos,
            // sparse: SparseTable::build(&block_mins),
            segtree: SegmentTree::build(block_mins.as_slice()),
        }
    }
    fn space(&self) -> usize {
        // self.sparse.space()
        self.segtree.space()
    }
    fn query(&self, data: impl Data<T>, l: usize, r: usize) -> (T, usize) {
        let block_l = l / S;
        let block_r = r / S;
        if block_l == block_r {
            return (l..r + 1).map(|i| (data.index(i), i)).min().unwrap();
        }
        let suffix = (l..(block_l + 1) * S)
            .map(|i| (data.index(i), i))
            .min()
            .unwrap();
        let prefix = (block_r * S..r + 1)
            .map(|i| (data.index(i), i))
            .min()
            .unwrap();
        let mid = if block_r > block_l + 1 {
            let (val, block_idx) = self.segtree.query(data, block_l + 1, block_r - 1);
            let idx = self.block_min_pos[block_idx];
            (val, block_idx * S + idx as usize)
        } else {
            (T::MAX, usize::MAX)
        };
        suffix.min(prefix).min(mid)
    }
}

/// Like `BlockRmq`, but precomputes prefix and suffix minima within each block so that
/// cross-block queries need only two O(1) lookups instead of two linear scans.
///
/// `prefix_min[i]` = min(data[block_start ..= i])
/// `suffix_min[i]` = min(data[i ..= block_end])
pub struct BlockRmqPrecomputed<T: RmqElem, const S: usize> {
    prefix_min: Vec<(T, usize)>,
    suffix_min: Vec<(T, usize)>,
    block_min_pos: Vec<usize>,
    // sparse: SparseTable<T>,
    segtree: SegmentTree<T>,
}
impl<T: RmqElem, const S: usize> Rmq<T> for BlockRmqPrecomputed<T, S> {
    fn name() -> String {
        format!("BlockPSRmq<{S}>")
    }
    fn build(data: impl Data<T>) -> Self {
        let n = data.lenn();
        let num_blocks = (n + S - 1) / S;
        let mut prefix_min = vec![(T::MAX, 0); n];
        let mut suffix_min = vec![(T::MAX, 0); n];

        for b in 0..num_blocks {
            let lo = b * S;
            let hi = (lo + S).min(n);
            // Prefix minima: left to right within block.
            prefix_min[lo] = (data.index(lo), lo);
            for i in lo + 1..hi {
                prefix_min[i] = prefix_min[i - 1].min((data.index(i), i));
            }
            // Suffix minima: right to left within block.
            suffix_min[hi - 1] = (data.index(hi - 1), hi - 1);
            for i in (lo..hi - 1).rev() {
                suffix_min[i] = suffix_min[i + 1].min((data.index(i), i));
            }
        }

        let (block_mins, block_min_pos): (Vec<T>, Vec<usize>) = (0..num_blocks)
            .into_par_iter()
            .map(|b| prefix_min[(b * S + S).min(n) - 1])
            .unzip();

        BlockRmqPrecomputed {
            prefix_min,
            suffix_min,
            block_min_pos,
            // sparse: SparseTable::build(&block_mins),
            segtree: SegmentTree::build(block_mins.as_slice()),
        }
    }
    fn space(&self) -> usize {
        std::mem::size_of_val(self.prefix_min.as_slice())
            + std::mem::size_of_val(self.suffix_min.as_slice())
            // + self.sparse.space()
            + self.segtree.space()
    }
    fn query(&self, data: impl Data<T>, l: usize, r: usize) -> (T, usize) {
        let block_l = l / S;
        let block_r = r / S;
        if block_l == block_r {
            return (l..=r).map(|i| (data.index(i), i)).min().unwrap();
        }
        let suffix = self.suffix_min[l];
        let prefix = self.prefix_min[r];
        let mid = if block_r > block_l + 1 {
            let (val, block_idx) = self.segtree.query(data, block_l + 1, block_r - 1);
            let idx = self.block_min_pos[block_idx];
            (val, idx)
        } else {
            (T::MAX, usize::MAX)
        };
        suffix.min(prefix).min(mid)
    }
}

type I = SaElem;

/// Sparse table: O(n log n) space, O(1) query.
///
/// `table[k][i]` = minimum of `data[i .. i + 2^k]`.
/// Query [l, r]: let k = floor(log2(r - l + 1)), return min(table[k][l], table[k][r - 2^k + 1]).
struct SparseTable<T: RmqElem> {
    table: Vec<Vec<(T, I)>>,
}
impl<T: RmqElem> Rmq<T> for SparseTable<T> {
    fn name() -> String {
        "SparseTable".to_string()
    }
    fn build(data: impl Data<T>) -> Self {
        let n = data.lenn();
        let levels = if n <= 1 {
            1
        } else {
            usize::BITS as usize - n.leading_zeros() as usize
        };

        let mut table: Vec<Vec<(T, I)>> = Vec::with_capacity(levels);
        table.push((0..n).map(|i| (data.index(i), i as I)).collect());
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
    fn query(&self, _data: impl Data<T>, l: usize, r: usize) -> (T, usize) {
        let k = usize::BITS as usize - (r - l + 1).leading_zeros() as usize - 1;
        let (x, idx) = self.table[k][l].min(self.table[k][r - (1 << k) + 1]);
        (x, idx as usize)
    }
}

/// Segment tree storing the range minimum.
///
/// Uses a 1-indexed implicit binary tree over the next power-of-two many leaves.
/// `tree[1]` is the root; children of node `i` are `2i` and `2i+1`.
/// Leaves are at indices `[size, 2*size)` where `size = n.next_power_of_two()`.
/// Space: 2 × (next power of two ≥ n) × 8 bytes ≈ 2n words.
/// Build: O(n).  Query: O(log n).
struct SegmentTree<T: RmqElem> {
    tree: Vec<(T, I)>,
    size: usize,
}

impl<'a, T: RmqElem> Rmq<T> for SegmentTree<T> {
    fn name() -> String {
        "SegmentTree".to_string()
    }
    fn build(data: impl Data<T>) -> Self {
        let n = data.lenn();
        let size = n.next_power_of_two();
        let mut tree = vec![(T::MAX, 0); 2 * size];
        // Copy leaves.
        for i in 0..n {
            tree[size + i] = (data.index(i), i as I);
        }
        // Fill internal nodes bottom-up.
        for i in (1..size).rev() {
            tree[i] = tree[2 * i].min(tree[2 * i + 1]);
        }
        SegmentTree { tree, size }
    }
    fn space(&self) -> usize {
        std::mem::size_of_val(self.tree.as_slice())
    }
    fn query(&self, _data: impl Data<T>, mut l: usize, mut r: usize) -> (T, usize) {
        let mut res = (T::MAX, I::MAX);
        l += self.size;
        r += self.size + 1;
        while l < r {
            if l & 1 == 1 {
                res = res.min(self.tree[l]);
                l += 1;
            }
            if r & 1 == 1 {
                r -= 1;
                res = res.min(self.tree[r]);
            }
            l >>= 1;
            r >>= 1;
        }
        let (x, idx) = res;
        (x, idx as usize)
    }
}
