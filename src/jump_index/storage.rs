#![cfg(feature = "mphf")]
use std::cmp::max;

use itertools::Itertools;
use mem_dbg::MemSize;
use ptr_hash::{hash::FastIntHash, DefaultPtrHash, PtrHash, PtrHashParams};
use sux::{bits::BitVec, rank_sel::SelectAdaptConst};

use crate::gbs;

use super::link::{BareEf, Link, LinkEf};

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
