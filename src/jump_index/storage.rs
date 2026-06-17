#![allow(unused_imports)]
use std::{cmp::max, collections::BTreeSet};

use itertools::{Either, Itertools};
use mem_dbg::MemSize;
#[cfg(feature = "mphf")]
use ptr_hash::{hash::FastIntHash, DefaultPtrHash, PtrHash, PtrHashParams};
use sux::{
    bits::BitVec,
    dict::EliasFanoBuilder,
    rank_sel::SelectAdaptConst,
    traits::{Pred, Succ},
};

use crate::gbs;

use super::link::{BareEf, Link, LinkEf, SuffixLink};

#[cfg(feature = "mphf")]
#[derive(MemSize)]
pub struct MphfStore {
    /// PHF mapping (source, c) to index.
    phf: DefaultPtrHash<FastIntHash, usize>,
    /// Bitvec indicating for each slot.
    /// 1 per slot
    /// 0 before the 1 for every additional link target here.
    bits: SelectAdaptConst<BitVec>,
    /// targets
    targets: Vec<u32>,
}

#[cfg(feature = "mphf")]
impl MphfStore {
    pub fn new(links: &LinkEf) -> Self {
        let keys: Vec<_> = links
            .iter()
            .map(|x| Link::from_key(x).source_c())
            .dedup()
            .collect();
        let mut params = PtrHashParams::default();
        params.remap = false;
        params.alpha = 0.98;
        eprintln!("Building PHF..");
        let phf = DefaultPtrHash::new(&keys, params);
        let slots = phf.max_index();
        let mut lens = vec![0u8; slots];
        for link in links.iter() {
            let link = Link::from_key(link);
            let idx = phf.index_no_remap(&link.source_c());
            lens[idx] += 1;
        }
        // increment all lens to at least 1
        for len in &mut lens {
            *len = 1.max(*len);
        }
        let target_slots = lens.iter().map(|x| *x as usize).sum();

        // build bitvec. 1 + 0*(len-1) for each slot, + 1 at the end.
        eprintln!("Build bitvec..");
        let mut bits = BitVec::with_capacity(lens.iter().map(|x| *x as usize).sum());
        for len in &lens {
            bits.push(true);
            bits.extend(std::iter::repeat(false).take((*len - 1) as usize));
        }
        bits.push(true);
        assert_eq!(bits.len(), target_slots + 1);
        let bits = SelectAdaptConst::<_, _>::new(bits);

        eprintln!("Build targets..");
        lens.fill(0);
        let mut targets = vec![0u32; target_slots];
        for link in links.iter() {
            let link = Link::from_key(link);
            let idx = phf.index_no_remap(&link.source_c());
            let pos = lens[idx] as usize;
            targets[pos] = link.target() as u32;
            lens[idx] += 1;
        }
        MphfStore { phf, bits, targets }
    }

    pub fn size(&self) -> String {
        format!(
            "MphfStore:\n  phf: {:.3} GB\n  bits: {:.3} GB\n  targets: {:.3} GB",
            gbs(&self.phf),
            gbs(&self.bits),
            gbs(&self.targets)
        )
    }
}

pub struct RelativeStore {
    fwd_ef: LinkEf,
    suf_ef: LinkEf,
    fwd_set: BTreeSet<u128>,
    suf_set: BTreeSet<u128>,
}

impl RelativeStore {
    pub fn new(fwd_ef: LinkEf, suf_ef: LinkEf) -> Self {
        Self {
            fwd_ef,
            suf_ef,
            fwd_set: BTreeSet::new(),
            suf_set: BTreeSet::new(),
        }
    }
    pub fn fwd_len(&self) -> usize {
        self.fwd_ef.len() + self.fwd_set.len()
    }
    pub fn suf_len(&self) -> usize {
        self.suf_ef.len() + self.suf_set.len()
    }
    pub fn insert_fwd(&mut self, link: Link) {
        self.fwd_set.insert(link.key());
    }
    pub fn insert_suf(&mut self, link: SuffixLink) {
        self.suf_set.insert(link.key());
    }
    /// Return the first fwd link >= `link`, but only if it has matching source and character.
    pub fn fwd_succ(&self, link: Link) -> Option<Link> {
        let x = link.key();
        let k1 = self.fwd_ef.succ(x).map(|x| x.1);
        let k2 = self.fwd_set.range(x..).next().copied();
        let k = match (k1, k2) {
            (Some(k1), Some(k2)) => Some(std::cmp::min(k1, k2)),
            (Some(k1), None) => Some(k1),
            (None, Some(k2)) => Some(k2),
            (None, None) => None,
        }?;
        let l = Link::from_key(k);
        if l.source_c() == link.source_c() {
            Some(l)
        } else {
            None
        }
    }
    /// Return the last fwd link <= `link`, but only if it has matching source and character.
    pub fn fwd_pred(&self, link: Link) -> Option<Link> {
        let x = link.key();
        let k1 = self.fwd_ef.pred(x).map(|x| x.1);
        let k2 = self.fwd_set.range(..=x).next_back().copied();
        let k = std::cmp::max(k1, k2)?;
        let l = Link::from_key(k);
        if l.source_c() == link.source_c() {
            Some(l)
        } else {
            None
        }
    }
    /// Return an iterator over suffix links >= link.
    pub fn suf_iter_from(&self, link: SuffixLink) -> impl Iterator<Item = SuffixLink> {
        let x = link.key();
        let k1 = self
            .suf_ef
            .iter_from_succ(x)
            .map(|x| Either::Left(x.1))
            .unwrap_or(Either::Right(std::iter::empty()));
        let k2 = self.suf_set.range(x..).copied();
        merging_iterator::MergeIter::new(k1, k2).map(SuffixLink::from_key)
    }
    pub fn finish(self) -> (LinkEf, LinkEf) {
        eprintln!("Merging fwd..");
        let fwd_ef = merge_ef_and_set(self.fwd_ef, self.fwd_set);
        eprintln!("Merging suf..");
        let suf_ef = merge_ef_and_set(self.suf_ef, self.suf_set);
        (fwd_ef, suf_ef)
    }
}

fn merge_ef_and_set(ef: LinkEf, set: BTreeSet<u128>) -> LinkEf {
    let n = ef.len() + set.len();
    let last = ef.upper_bound().max(set.last().copied().unwrap_or(0));
    let mut builder = EliasFanoBuilder::new(n, last);
    for k in merging_iterator::MergeIter::new(ef.into_iter(), set.into_iter()) {
        builder.push(k);
    }
    let fwd_ef = builder.build_with_dict();
    fwd_ef
}
