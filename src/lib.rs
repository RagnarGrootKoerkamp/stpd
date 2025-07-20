#![feature(gen_blocks)]
use std::{cmp::{Ordering, Reverse}, collections::{HashMap, HashSet}};
use itertools::Itertools;
use rand::{rng, seq::SliceRandom};

pub mod strings;
pub mod test;

pub type T = Vec<u8>;
pub type SA = Vec<usize>;
pub type LCP = Vec<usize>;

pub fn print(t: &[u8]) -> String {
    String::from_utf8(t.iter().map(|c| b'0' + c).collect_vec()).unwrap()
}

pub fn sa(t: &T) -> SA {
    let mut sa = (0..t.len()).collect_vec();
    sa.sort_by_key(|&i| &t[i..]);
    sa
}

fn co_sa(t: &T) -> Vec<usize> {
    let mut co_sa = (0..t.len()).collect_vec();
    co_sa.sort_by(|&i, &j|{
        let mut i = i + 1;
        let mut j = j + 1;
        while i > 0 && j > 0{
            if t[i-1] != t[j-1] {
                return t[i-1].cmp(&t[j-1]);
            }
            i -= 1;
            j -= 1;
        }
        if i == 0 {
            return Ordering::Less;
        }
        if j == 0 {
            return Ordering::Greater;
        }
        unreachable!()
    });
    co_sa
}

pub fn lcp(t: &T, sa: &SA) -> LCP {
    let n = t.len();
    let lcp = |(i, j)| {
        let mut l = 0;
        while i + l < n && j + l < n && t[i + l] == t[j + l] {
            l += 1;
        }
        l
    };
    sa.iter().circular_tuple_windows().map(lcp).collect()
}

pub fn bwt(t: &T, sa: &SA) -> T {
    sa.iter()
        .map(|&i| t[if i == 0 { t.len() - 1 } else { i - 1 }])
        .collect_vec()
}

/// Number of BWT runs.
pub fn r(bwt: &T) -> usize {
    1 + bwt
        .iter()
        .circular_tuple_windows()
        .map(|(l, r)| if l != r { 1 } else { 0 })
        .sum::<usize>()
}

/// delta, and maximizing k.
pub fn delta(t: &T) -> (f32, usize) {
    let mut max = (0.0, 0);
    for k in 1..=t.len() {
        let mut kmers = HashSet::new();
        for kmer in t.windows(k) {
            kmers.insert(kmer);
        }
        let delta_k = (kmers.len() as f32 / k as f32, k);
        if delta_k > max {
            max = delta_k;
        }
    }
    max
}

/// Iterate right maximal strings α, corresponding to paths to explicit suffix tree nodes.
pub fn tree_nodes<'t>(t: &'t T, sa: &SA, lcp: &LCP) -> impl Iterator<Item = &'t [u8]> {
    gen {
        let mut depths = vec![0];
        for i in 0..t.len() {
            let l = lcp[i];
            while *depths.last().unwrap() > l {
                let l2 = depths.pop().unwrap();
                yield &t[sa[i]..sa[i]+l2];
            }
            if *depths.last().unwrap() < l {
                depths.push(l);
            }
        }
        yield b"";
    }
}

/// Iterate right maximal extensions αX, corresponding to edges going out of explicit suffix tree nodes.
pub fn tree_edges<'t>(t: &'t T, sa: &SA, lcp: &LCP) -> impl Iterator<Item = &'t [u8]> {
    gen {
        let mut depths = vec![];
        for i in 0..t.len() {
            let l = lcp[i];
            // eprintln!("{i} {l} {} {:?}", sa[i], &t[sa[i]..]);
            while let Some(l2) = depths.last() && *l2 > l {
                let l2 = depths.pop().unwrap();
                yield &t[sa[i]..sa[i]+l2+1];
            }
            // Otherwise the suffix is just an internal node anyway.
            if sa[i] + l + 1 <=  t.len(){
                yield &t[sa[i]..sa[i]+l+1];
            }
            depths.push(l);
        }
    }
}

/// Number of leaves in the Weiner-link-tree on explicit suffix tree nodes.
pub fn w(t: &T, sa: &SA, lcp: &LCP) -> usize {
    let mut internal = HashSet::new();
    let mut num_nodes = 0;
    for n in tree_nodes(t,sa,lcp) {
        num_nodes += 1;
        if !n.is_empty() {
            internal.insert(&n[1..]);
        }
    }
    num_nodes - internal.len()
}

/// Size of smallest suffixient set.
pub fn chi(t: &T, sa: &SA, lcp: &LCP, print: bool) -> usize {
    let mut edges = tree_edges(t, sa, lcp).collect_vec();
    edges.sort_by_key(|e| Reverse(e.len()));

    let mut chi = 0;
    let mut covered = HashSet::new();
    covered.insert([].as_slice());
    for mut e in edges {
        if !covered.contains(e) {
            chi += 1;

            if print {
                eprintln!("{}", crate::print(e));
            }

            while !covered.contains(e) {
                covered.insert(e);
                e = &e[1..];
            }
        }
    }
    chi
}

// Returns a map alphaX -> alphaX..., where each right-maximal extension is mapped to its node.
pub fn tree<'t>(t: &'t T, sa: &SA, lcp: &LCP) -> HashMap<&'t [u8], &'t[u8]> {
    let mut nodes = HashSet::new();
    for n in tree_nodes(t, sa, lcp) {
        nodes.insert(n);
    }
    let mut tree = HashMap::new();
    for e in tree_edges(t, sa, lcp) {
        // Iterate up until a node is found.
        let mut p = e;
        while !nodes.contains(p.split_last().unwrap().1) {
            p.split_off_last().unwrap();
        }
        tree.insert(p, e);
    }
    tree
}

/// Size of smallest suffixient set, with path decomposition.
pub fn chi_pd(t: &T, sa: &SA, lcp: &LCP) -> usize {
    let mut edges = tree_edges(t, sa, lcp).collect_vec();
    edges.sort_by_key(|e| Reverse(e.len()));

    let tree = tree(t, sa, lcp);

    let mut chi = 0;
    let mut covered = HashSet::new();
    let mut path_covered = HashSet::new();
    covered.insert([].as_slice());
    for mut e in edges {

        if !covered.contains(e) && !path_covered.contains(e) {
            chi += 1;

            while !covered.contains(e) {
                covered.insert(e);

                // Also cover right-extensions of e, after walking down the tree.
                {
                    // Position of start of e.
                    let pos = unsafe { e.as_ptr().offset_from(t.as_ptr()) as usize };
                    let mut e = &e[..e.len()-1];
                    while pos + e.len()+1 < t.len() {
                        // right-extend by 1 char to get to a right-maximal-extension.
                        e = &t[pos..pos+e.len()+1];
                        path_covered.insert(e);
                        if let Some(e2) = tree.get(e) {
                            e = e2;
                        } else {
                            break;
                        }

                    }
                }


                e = &e[1..];
            }
        }
    }
    chi
}

/// Size of smallest suffixient set, with path decomposition.
pub fn chi_pd2(t: &T, sa: &SA, lcp: &LCP) -> usize {
    let mut edges = tree_edges(t, sa, lcp).collect_vec();
    edges.sort_by_key(|e| Reverse(e.len()));

    let tree = tree(t, sa, lcp);

    let mut chi = 0;
    let mut covered = HashSet::new();
    let mut path_covered = HashSet::new();
    covered.insert([].as_slice());
    for mut e in edges {

        if !covered.contains(e) && !path_covered.contains(e) {
            chi += 1;

            while e != b"" {
                covered.insert(e);

                // Also cover right-extensions of e, after walking down the tree.
                {
                    // Position of start of e.
                    let pos = unsafe { e.as_ptr().offset_from(t.as_ptr()) as usize };
                    let mut e = &e[..e.len()-1];
                    while pos + e.len()+1 < t.len() {
                        // right-extend by 1 char to get to a right-maximal-extension.
                        e = &t[pos..pos+e.len()+1];
                        path_covered.insert(e);
                        if let Some(e2) = tree.get(e) {
                            e = e2;
                        } else {
                            break;
                        }

                    }
                }


                e = &e[1..];
            }
        }
    }
    chi
}

pub fn stpd(t: &T, _sa: &SA, _lcp: &LCP, perm: &Vec<usize>) -> usize {
    let mut iperm = vec![0; t.len()];
    for i in 0..t.len(){
        iperm[perm[i]] = i;
    }
    let mut seen = HashSet::new();
    let mut sampled = HashSet::new();

    for len in 1..=t.len(){
        for &pos in &iperm {
            if pos < len-1 {continue;}
            let start = pos-len+1;
            let e = &t[start..=pos];
            if seen.contains(e) {
                continue;
            }
            sampled.insert(pos);
            for end in pos+1..=t.len(){
                seen.insert(&t[start..end]);
            }
        }
    }
    sampled.len()
}

pub fn stpd_pos_minus(t: &T, sa: &SA, lcp: &LCP) -> usize {
    let perm = (0..t.len()).collect_vec();
    stpd(t, sa, lcp, &perm)
}

pub fn stpd_pos_plus(t: &T, sa: &SA, lcp: &LCP) -> usize {
    let perm = (0..t.len()).rev().collect_vec();
    stpd(t, sa, lcp, &perm)
}

pub fn stpd_lex_minus(t: &T, sa: &SA, lcp: &LCP) -> usize {
    let mut isa = vec![0; t.len()];
    for (i, &x) in sa.iter().enumerate(){
        isa[x] = i;
    }
    stpd(t, sa, lcp, &isa)
}

pub fn stpd_lex_plus(t: &T, sa: &SA, lcp: &LCP) -> usize {
    let mut isa = vec![0; t.len()];
    for (i, &x) in sa.iter().enumerate(){
        isa[x] = t.len()-1-i;
    }
    stpd(t, sa, lcp, &isa)
}

pub fn stpd_colex_minus(t: &T, sa: &SA, lcp: &LCP) -> usize {
    let co_sa = co_sa(t);
    let mut i_co_sa = vec![0; t.len()];
    for (i, &x) in co_sa.iter().enumerate(){
        i_co_sa[x] = i;
    }
    stpd(t, sa, lcp, &i_co_sa)
}

pub fn stpd_colex_plus(t: &T, sa: &SA, lcp: &LCP) -> usize {
    let co_sa = co_sa(t);
    let mut i_co_sa = vec![0; t.len()];
    for (i, &x) in co_sa.iter().enumerate(){
        i_co_sa[x] = t.len()-1-i;
    }
    stpd(t, sa, lcp, &i_co_sa)
}


pub fn stpd_rand(t: &T, sa: &SA, lcp: &LCP) -> usize {
    let mut perm = (0..t.len()).collect_vec();
    perm.shuffle(&mut rng());
    stpd(t, sa, lcp, &perm)
}
