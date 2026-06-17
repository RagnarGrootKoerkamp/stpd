#![allow(unreachable_code, unused)]
use std::cmp::Ordering::{Greater, Less};

use crate::{
    rmq::{self, Rmq},
    stpd::cmp_colex,
    T,
};

use super::{link, JumpIndex, Pi};

pub struct Stpd<'t, const PI: Pi> {
    pub t: &'t [u8],
    pub stpd_samples: Vec<usize>,
    pub stpd_pi: Vec<u64>,
    pub stpd_rmq: rmq::BlockRmq<u64, 64>,
}

impl<'t, const PI: Pi> Stpd<'t, PI> {
    pub fn from_jump_index(ji: &JumpIndex<'t, PI>) -> Self {
        let mut stpd_samples: Vec<usize> = ji
            .fwd_links
            .iter()
            .map(|l| link::Link::from_key(l).target())
            .collect();
        stpd_samples.sort_by(|&a, &b| cmp_colex(&ji.t[..=a], &ji.t[..=b]).1);
        stpd_samples.dedup();
        let _stpd_pi = todo!("THIS NEEDS THE SUFFIX ARRAY FOR Pi::LeftMost.");
        // let stpd_pi: Vec<u64> = stpd_samples.iter().map(|&x| pi[x] as u64).collect();
        Stpd {
            t: ji.t,
            stpd_samples,
            stpd_pi: _stpd_pi,
            stpd_rmq: rmq::BlockRmq::build(_stpd_pi.as_slice()),
        }
    }

    /// Returns leftmost text position where the pattern matches.
    pub fn map_stpd(&self, pattern: &[u8]) -> Option<usize> {
        let mut pos = 0;
        for (i, &c) in pattern.iter().enumerate() {
            if self.t[pos] == c {
                pos += 1;
                continue;
            }

            // TODO: Use binary search function that reuses LCP.
            let idx1 = self
                .stpd_samples
                .binary_search_by(|&sample_pos| {
                    match cmp_colex(&self.t[..=sample_pos], &pattern[..=i]).1 {
                        std::cmp::Ordering::Equal => Greater,
                        x => x,
                    }
                })
                .unwrap_err();
            let idx2 = self
                .stpd_samples
                .binary_search_by(|&sample_pos| {
                    match cmp_colex(&self.t[..=sample_pos], &pattern[..=i]).1 {
                        std::cmp::Ordering::Equal => Less,
                        x => x,
                    }
                })
                .unwrap_err();
            if idx1 == idx2 {
                // eprintln!("idx: {idx1}..{idx2}");
                return None;
            }
            let (_val, idx) = self.stpd_rmq.query(self.stpd_pi.as_slice(), idx1, idx2 - 1);
            // eprintln!("idx: {idx1}..{idx2} => {idx} val={_val}");
            pos = self.stpd_samples[idx] + 1;

            // eprintln!("pos: {pos}");
        }
        // eprintln!("end pos {pos}");
        Some(pos - pattern.len())
    }

    pub fn space(&self) {
        eprintln!(
            "stpd samples {:.3} GB",
            std::mem::size_of_val(self.stpd_samples.as_slice()) as f32 / 1e9
        );
        eprintln!(
            "stpd pi      {:.3} GB",
            std::mem::size_of_val(self.stpd_pi.as_slice()) as f32 / 1e9
        );
        eprintln!("stpd rmq     {:.3} GB", self.stpd_rmq.space() as f32 / 1e9);
    }
}
