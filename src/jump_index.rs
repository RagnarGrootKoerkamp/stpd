use itertools::Itertools;
use std::{
    cmp::Ordering::{self, Greater, Less},
    collections::BTreeSet,
    marker::Sync,
};

use crate::{
    bwt, print,
    rmq::{self, Rmq},
    sa_and_lcp,
    stpd::cmp_colex,
    LCP, SA, T,
};

#[derive(Copy, Clone, PartialEq, PartialOrd, Eq, Ord, Debug)]
pub struct Link {
    pub source: usize,
    pub c: u8,
    pub lcp: usize,
    pub target: usize,
}

impl Link {
    fn key(&self) -> u128 {
        ((self.source as u128) << 64) | ((self.c as u128) << 48) | (self.lcp as u128)
    }
}

impl voracious_radix_sort::Radixable<u128> for Link {
    type Key = u128;
    fn key(&self) -> Self::Key {
        self.key()
    }
}

pub struct JumpIndex<TR: AsRef<T>, SAR: AsRef<SA>, LCPR: AsRef<LCP>> {
    pub t: TR,
    pub sa: SAR,
    pub lcp: LCPR,
    pub stpd_samples: Vec<usize>,
    pub stpd_pi: Vec<u64>,
    pub stpd_rmq: rmq::BlockRmq<64>,
    // TODO: Predecessor structure
    pub links: Vec<Link>,
}

pub struct JumpIndexStats {
    pub num_sampled: usize,
    pub num_sources: usize,
    pub num_source_chars: usize,
    pub num_links: usize,
}

impl<TR: AsRef<T> + Sync> JumpIndex<TR, Vec<usize>, Vec<usize>> {
    pub fn new(t: TR) -> Self {
        let (sa, lcp) = sa_and_lcp(t.as_ref());
        let bwt = &bwt(t.as_ref(), &sa);
        let pi = (0..t.as_ref().len()).collect_vec();
        Self::new2(t, sa, bwt, lcp, &pi)
    }
}
impl<TR: AsRef<T> + Sync, SAR: AsRef<SA> + Sync, LCPR: AsRef<LCP> + Sync> JumpIndex<TR, SAR, LCPR> {
    pub fn new2(t: TR, sa: SAR, bwt: &T, lcp: LCPR, pi: &Vec<usize>) -> Self {
        const PARALLEL_THRESHOLD: usize = 100_000;

        struct State<'a, T, SA> {
            t: T,
            sa: SA,
            lcp: Vec<u64>,
            run_boundaries: BTreeSet<usize>,
            lcp_rmq: rmq::BlockRmq<128>,
            #[allow(unused)]
            permuted_pi: &'a Vec<u64>,
            pi_rmq: rmq::BlockRmq<128>,
        }

        impl<'a, TR: AsRef<T>, SAR: AsRef<SA>> State<'a, TR, SAR> {
            fn split(
                &self,
                interval: std::ops::Range<usize>,
            ) -> Option<(usize, usize, Vec<std::ops::Range<usize>>)> {
                if interval.len() <= 1 {
                    return None;
                }
                // FIXME
                // if self.run_boundaries.range(interval.clone()).next().is_none() {
                //     return None;
                // }

                let anchor_pos = self
                    .pi_rmq
                    .query(&self.permuted_pi, interval.start, interval.end - 1)
                    .1;
                let mut done_intervals = vec![];
                let mut wip_intervals = vec![interval.clone()];
                let lcp = self.lcp[self
                    .lcp_rmq
                    .query(&self.lcp, interval.start, interval.end - 2)
                    .1];
                while let Some(interval) = wip_intervals.pop() {
                    if interval.len() <= 1 {
                        done_intervals.push(interval);
                        continue;
                    }
                    let split_pos = self
                        .lcp_rmq
                        .query(&self.lcp, interval.start, interval.end - 2)
                        .1
                        + 1;
                    let new_lcp = self.lcp[split_pos - 1];
                    if new_lcp > lcp {
                        done_intervals.push(interval);
                        continue;
                    }
                    assert!(new_lcp == lcp);
                    wip_intervals.push(interval.start..split_pos);
                    wip_intervals.push(split_pos..interval.end);
                }
                Some((anchor_pos, lcp as usize, done_intervals))
            }

            fn node_output(
                &self,
                anchor_pos: usize,
                lcp: usize,
                done_intervals: &[std::ops::Range<usize>],
                sampled: &mut Vec<usize>,
                links: &mut Vec<Link>,
            ) {
                for x in done_intervals {
                    if !x.contains(&anchor_pos) {
                        let secondary_anchor_pos =
                            self.pi_rmq.query(&self.permuted_pi, x.start, x.end - 1).1;
                        let text_idx = self.sa.as_ref()[secondary_anchor_pos];
                        let target = text_idx + lcp;
                        if target < self.t.as_ref().len() {
                            sampled.push(target);
                            let source = self.sa.as_ref()[anchor_pos] + lcp;
                            let c = self.t.as_ref()[target];
                            links.push(Link {
                                source,
                                c,
                                lcp: co_lcp(&self.t.as_ref()[..source], &self.t.as_ref()[..target]),
                                target,
                            });
                        }
                    }
                }
            }

            fn dfs(
                &self,
                interval: std::ops::Range<usize>,
                sampled: &mut Vec<usize>,
                links: &mut Vec<Link>,
            ) {
                let Some((anchor_pos, lcp, done_intervals)) = self.split(interval) else {
                    return;
                };
                self.node_output(anchor_pos, lcp, &done_intervals, sampled, links);
                for x in done_intervals {
                    self.dfs(x, sampled, links);
                }
            }

            fn collect_work(
                &self,
                interval: std::ops::Range<usize>,
                sampled: &mut Vec<usize>,
                links: &mut Vec<Link>,
                work_queue: &mut Vec<std::ops::Range<usize>>,
            ) {
                if interval.len() < PARALLEL_THRESHOLD {
                    work_queue.push(interval);
                    return;
                }
                let Some((anchor_pos, lcp, done_intervals)) = self.split(interval) else {
                    return;
                };
                self.node_output(anchor_pos, lcp, &done_intervals, sampled, links);
                for x in done_intervals {
                    self.collect_work(x, sampled, links, work_queue);
                }
            }
        }

        let run_boundaries = (0..t.as_ref().len() - 1)
            .tuple_windows()
            .filter(|(i, j)| bwt[*i] != bwt[*j])
            .map(|(i, _j)| i)
            .collect_vec();
        let run_boundaries = BTreeSet::from_iter(run_boundaries);

        let permuted_pi: Vec<usize> = sa.as_ref().par_iter().map(|&i| pi[i]).collect();

        use rmq::Rmq as _;
        let lcp_u64: Vec<u64> = lcp.as_ref().iter().map(|&x| x as u64).collect();
        let permuted_pi: Vec<u64> = permuted_pi.iter().map(|&x| x as u64).collect();
        let state = State {
            t,
            sa,
            run_boundaries,
            lcp_rmq: rmq::BlockRmq::build(&lcp_u64),
            lcp: lcp_u64,
            pi_rmq: rmq::BlockRmq::build(&permuted_pi),
            permuted_pi: &permuted_pi,
        };

        let mut stpd_samples = vec![];
        let mut links = vec![];
        let mut work_queue = vec![];
        state.collect_work(
            0..state.t.as_ref().len(),
            &mut stpd_samples,
            &mut links,
            &mut work_queue,
        );

        use rayon::prelude::*;
        let child_results: Vec<(Vec<usize>, Vec<Link>)> = work_queue
            .into_par_iter()
            .map(|interval| {
                let mut sampled = vec![];
                let mut links = vec![];
                state.dfs(interval, &mut sampled, &mut links);
                (sampled, links)
            })
            .collect();
        for (s, l) in child_results {
            stpd_samples.extend(s);
            links.extend(l);
        }

        let State { t, sa, .. } = state;

        use voracious_radix_sort::RadixSort;
        stpd_samples.voracious_mt_sort(12);
        links.voracious_mt_sort(12);

        stpd_samples.dedup();
        links.dedup();

        stpd_samples.sort_by(|&a, &b| cmp_colex(&t.as_ref()[..=a], &t.as_ref()[..=b]).1);
        let stpd_pi: Vec<u64> = stpd_samples.iter().map(|&x| pi[x] as u64).collect();
        stpd_samples.iter().take(10).for_each(|&i| {
            eprintln!(
                "{i:>3}: {} ({})",
                print(&t.as_ref()[i.saturating_sub(30)..=i]),
                pi[i]
            );
        });

        JumpIndex {
            t,
            sa,
            lcp,
            stpd_samples,
            stpd_rmq: rmq::BlockRmq::build(&stpd_pi),
            stpd_pi,
            links,
        }
    }

    pub fn stats(&self) -> JumpIndexStats {
        JumpIndexStats {
            num_sampled: 1 + self
                .stpd_samples
                .iter()
                .tuple_windows()
                .filter(|(a, b)| a != b)
                .count(),
            num_sources: 1 + self
                .links
                .iter()
                .tuple_windows()
                .filter(|(a, b)| a.source != b.source)
                .count(),
            num_source_chars: 1 + self
                .links
                .iter()
                .tuple_windows()
                .filter(|(a, b)| (a.source, a.c) != (b.source, b.c))
                .count(),
            num_links: 1 + self
                .links
                .iter()
                .tuple_windows()
                .filter(|(a, b)| a != b)
                .count(),
        }
    }

    /// Returns leftmost text position where the pattern matches.
    pub fn map_stpd(&self, pattern: &[u8]) -> Option<usize> {
        let mut pos = 0;
        for (i, &c) in pattern.iter().enumerate() {
            if self.t.as_ref()[pos] == c {
                pos += 1;
                continue;
            }

            // TODO: Use binary search function that reuses LCP.
            let idx1 = self
                .stpd_samples
                .binary_search_by(|&sample_pos| {
                    match cmp_colex(&self.t.as_ref()[..=sample_pos], &pattern[..=i]).1 {
                        std::cmp::Ordering::Equal => Greater,
                        x => x,
                    }
                })
                .unwrap_err();
            let idx2 = self
                .stpd_samples
                .binary_search_by(|&sample_pos| {
                    match cmp_colex(&self.t.as_ref()[..=sample_pos], &pattern[..=i]).1 {
                        std::cmp::Ordering::Equal => Less,
                        x => x,
                    }
                })
                .unwrap_err();
            if idx1 == idx2 {
                // eprintln!("idx: {idx1}..{idx2}");
                return None;
            }
            let (_val, idx) = self.stpd_rmq.query(&self.stpd_pi, idx1, idx2 - 1);
            // eprintln!("idx: {idx1}..{idx2} => {idx} val={_val}");
            pos = self.stpd_samples[idx] + 1;

            // eprintln!("pos: {pos}");
        }
        // eprintln!("end pos {pos}");
        Some(pos - pattern.len())
    }

    /// Returns leftmost text position where the pattern matches.
    pub fn map_jump(&self, pattern: &[u8]) -> Option<usize> {
        let mut pos = 0;
        for (i, &c) in pattern.iter().enumerate() {
            if self.t.as_ref()[pos] == c {
                pos += 1;
                continue;
            }

            let link_idx = self
                .links
                .binary_search_by(|link| {
                    link.key().cmp(
                        &Link {
                            source: pos,
                            c,
                            lcp: i,
                            target: 0,
                        }
                        .key(),
                    )
                })
                .map_or_else(|e| e, |v| v);
            let link = self.links[link_idx];
            // eprintln!("pos {pos} link {link:?}");
            if link.source == pos && link.c == c {
                pos = link.target + 1;
            } else {
                return None;
            }
        }
        Some(pos - pattern.len())
    }
}

fn co_lcp(a: &[u8], b: &[u8]) -> usize {
    let min = a.len().min(b.len());
    for i in 0..min {
        if a[a.len() - 1 - i] != b[b.len() - 1 - i] {
            return i;
        }
    }
    min
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::strings::relative;

    #[test]
    fn test() {
        let mut t1 = std::time::Duration::default();
        let mut t2 = std::time::Duration::default();
        for (len, repeats, r) in [
            (10, 1, 0.1),
            (10, 2, 0.1),
            (10, 3, 0.1),
            (10, 4, 0.1),
            (10, 10, 0.1),
            (100, 1, 0.1),
            (100, 10, 0.1),
            (100, 100, 0.05),
            (100, 1000, 0.05),
            (1000, 100, 0.05),
            (10000, 10, 0.05),
        ] {
            let t = relative(len, 4, repeats, r).1;
            eprintln!("text: {}", print(&t));

            eprintln!("building for {len}x{repeats} at {r}..");
            let ji = JumpIndex::new(&t);

            // find a bunch of random substrings
            eprintln!("querying..");
            for len in 0..=len.min(1000) {
                let pos = rand::random_range(0..=t.len() - len);
                let pattern = &t[pos..pos + len];

                eprintln!("pattern: {}", print(pattern));

                let s = std::time::Instant::now();
                let p1 = ji.map_stpd(pattern);
                t1 += s.elapsed();
                let s = std::time::Instant::now();
                let p2 = ji.map_jump(pattern);
                t2 += s.elapsed();
                eprintln!("p1: {p1:?}");
                eprintln!("p2: {p2:?}");
                let p1 = p1.unwrap();
                let p2 = p2.unwrap();

                assert_eq!(&t[p1..p1 + len], pattern);
                assert_eq!(&t[p2..p2 + len], pattern);
                assert_eq!(p1, p2);
            }
        }
        eprintln!("STPD: {t1:?}");
        eprintln!("JI:   {t2:?}");
    }
}
