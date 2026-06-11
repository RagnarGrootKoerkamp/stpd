use itertools::Itertools;
use sux::bits::BitVec;

use sux::rank_sel::SelectZeroAdaptConst;

use sux::dict::{EliasFano, EliasFanoBuilder};
use voracious_radix_sort::RadixSort;

#[derive(Copy, Clone, PartialEq, Eq)]
pub struct Link {
    pub data: u128,
    // source: usize,
    // c: u8,
    // lcp: usize,
    // target: usize,
}

/// Variant without LCP value, as most (src, c) only have 1 target anyway.
#[derive(Copy, Clone, PartialEq, Eq)]
pub struct CompactLink {
    pub data: u128,
    // source: usize,
    // c: u8,
    // target: usize,
}

pub const SOURCE_BITS: u32 = 32;
pub const C_BITS: u32 = 8; // TODO: Reduce to 2
pub const LCP_BITS: u32 = 22; // enough for reference genome
pub const TARGET_BITS: u32 = 32; // enough for reference genome
pub const LINK_BITS: u32 = SOURCE_BITS + C_BITS + LCP_BITS + TARGET_BITS;

pub type LinkEf = EliasFano<u128, SelectZeroAdaptConst<BitVec<Box<[usize]>>, Box<[usize]>, 12, 3>>;
pub type BareEf = EliasFano<u128>;

impl Link {
    // const MAX: u128 = (1 << LINK_BITS) - 1;
    pub fn from_key(data: u128) -> Self {
        Self { data }
    }
    pub fn new(source: usize, c: u8, lcp: usize, target: usize) -> Self {
        assert!(LINK_BITS <= 128);
        assert!(
            source < (1 << SOURCE_BITS),
            "link {source},{c} -> {lcp},{target}"
        );
        assert!(
            (c as usize) < (1 << C_BITS),
            "link {source},{c} -> {lcp},{target}"
        );
        assert!(lcp < (1 << LCP_BITS), "link {source},{c} -> {lcp},{target}");
        assert!(
            target < (1 << TARGET_BITS),
            "link {source},{c} -> {lcp},{target}"
        );

        let data = ((source as u128) << (C_BITS + LCP_BITS + TARGET_BITS))
            | ((c as u128) << (LCP_BITS + TARGET_BITS))
            | ((lcp as u128) << TARGET_BITS)
            | (target as u128);
        Self { data }
    }
    pub fn key(&self) -> u128 {
        self.data
    }
    pub fn source(&self) -> usize {
        ((self.key() >> (C_BITS + LCP_BITS + TARGET_BITS)) & ((1 << SOURCE_BITS) - 1)) as usize
    }
    pub fn source_c(&self) -> usize {
        ((self.key() >> (LCP_BITS + TARGET_BITS)) & ((1 << (SOURCE_BITS + C_BITS)) - 1)) as usize
    }
    pub fn c(&self) -> u8 {
        ((self.key() >> (LCP_BITS + TARGET_BITS)) & ((1 << C_BITS) - 1)) as u8
    }
    pub fn lcp(&self) -> usize {
        ((self.key() >> TARGET_BITS) & ((1 << LCP_BITS) - 1)) as usize
    }
    pub fn target(&self) -> usize {
        (self.key() & ((1 << TARGET_BITS) - 1)) as usize
    }

    pub fn links_to_ef(links: Vec<Link>) -> BareEf {
        let mut links = links.into_iter().map(|l| l.key()).collect_vec();
        links.voracious_sort();
        links.dedup();
        BareEf::from(links)
    }

    pub fn compactify(&self) -> CompactLink {
        CompactLink::new(self.source(), self.c(), self.target())
    }

    /// One EF with (src, c, target) with a single target per (src, c),
    /// and a second one with all (src, c, lcp, target) (src, c) with additional targets.
    pub fn links_to_compact_ef(links: Vec<Link>) -> (BareEf, BareEf) {
        let num_sources = links.chunk_by(|x, y| x.source_c() == y.source_c()).count();
        let mut links = links.into_iter().map(|l| l.key()).collect_vec();
        links.voracious_sort();
        links.dedup();

        let u_full = *links.last().unwrap();
        let u_compact = Link::from_key(u_full).compactify().key();
        let mut ef_compact = EliasFanoBuilder::new(num_sources, u_compact);
        let mut ef_full = EliasFanoBuilder::new(links.len() - num_sources, u_full);

        for chunk in
            links.chunk_by(|x, y| Link::from_key(*x).source_c() == Link::from_key(*y).source_c())
        {
            ef_compact.push(Link::from_key(chunk[0]).compactify().key());
            for &link in &chunk[1..] {
                ef_full.push(link);
            }
        }
        (ef_compact.build(), ef_full.build())
    }
}

impl std::fmt::Debug for Link {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Link")
            .field("source", &self.source())
            .field("c", &self.c())
            .field("lcp", &self.lcp())
            .field("target", &self.target())
            .finish()
    }
}

impl CompactLink {
    // const MAX: u128 = (1 << LINK_BITS) - 1;
    pub fn from_key(data: u128) -> Self {
        Self { data }
    }
    pub fn new(source: usize, c: u8, target: usize) -> Self {
        assert!(LINK_BITS <= 128);
        assert!(source < (1 << SOURCE_BITS), "link {source},{c} -> {target}");
        assert!(
            (c as usize) < (1 << C_BITS),
            "link {source},{c} -> {target}"
        );
        assert!(target < (1 << TARGET_BITS), "link {source},{c} -> {target}");
        let data = ((source as u128) << (C_BITS + TARGET_BITS))
            | ((c as u128) << (TARGET_BITS))
            | (target as u128);
        Self { data }
    }
    pub fn key(&self) -> u128 {
        self.data
    }
    pub fn source(&self) -> usize {
        ((self.key() >> (C_BITS + TARGET_BITS)) & ((1 << SOURCE_BITS) - 1)) as usize
    }
    pub fn source_c(&self) -> usize {
        ((self.key() >> TARGET_BITS) & ((1 << (SOURCE_BITS + C_BITS)) - 1)) as usize
    }
    pub fn c(&self) -> u8 {
        ((self.key() >> TARGET_BITS) & ((1 << C_BITS) - 1)) as u8
    }
    pub fn target(&self) -> usize {
        (self.key() & ((1 << TARGET_BITS) - 1)) as usize
    }

    /// Store all values 'mirrored': `u-x` for x from large to small,
    /// where `u` is the largest element.
    pub fn links_to_ef(links: Vec<Link>) -> BareEf {
        let mut links = links.into_iter().map(|l| l.key()).collect_vec();
        links.voracious_sort();
        links.dedup();
        BareEf::from(links)
    }
}

impl std::fmt::Debug for CompactLink {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Link")
            .field("source", &self.source())
            .field("c", &self.c())
            .field("target", &self.target())
            .finish()
    }
}
