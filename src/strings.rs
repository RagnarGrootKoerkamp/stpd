use crate::T;
use itertools::Itertools;
use packed_seq::{AsciiSeqVec, Seq, SeqVec};
use rand::seq::IteratorRandom;

pub fn random(n: usize, sigma: u8) -> (String, T) {
    let mut t = (0..n).map(|_| rand::random_range(1..=sigma)).collect_vec();
    t.push(0);
    (format!("rand({n},{sigma})"), t)
}

pub fn relative(n: usize, sigma: u8, copies: usize, r: f32) -> (String, T) {
    let reference = (0..n).map(|_| rand::random_range(1..=sigma)).collect_vec();
    let mutations = (n as f32 * r).ceil() as usize;
    let mut t = vec![];
    for _ in 0..copies {
        let mut copy = reference.clone();
        for i in (0..n).choose_multiple(&mut rand::rng(), mutations) {
            let mut c = rand::random_range(1..=sigma - 1);
            if c >= copy[i] {
                c += 1;
            }
            copy[i] = c;
        }
        t.extend(copy);
    }
    t.push(0);
    (format!("relative({copies}*{n}@{r},{sigma})"), t)
}

pub fn fib(n: usize) -> (String, T) {
    let mut t = vec![1, 2];
    let mut last = 1;
    let mut cur = 2;
    for _ in 2..n {
        t.extend_from_within(..last);

        let new = last + cur;
        last = cur;
        cur = new;
    }
    (format!("fib({n})"), t)
}

pub fn thue_morse(n: usize) -> (String, T) {
    let mut t = vec![0; 1 << n];
    t[0] = 1;

    for i in 0..n {
        let l = 1 << i;
        for j in 0..l {
            t[l + j] = 3 - t[j];
        }
    }
    t.push(0);

    (format!("thue_morse({n})"), t)
}

pub fn rev((n, mut t): (String, T)) -> (String, T) {
    t.reverse();
    (format!("rev@{n}"), t)
}

pub fn flip((n, mut t): (String, T)) -> (String, T) {
    for x in &mut t {
        *x = 3 - *x;
    }
    (format!("flip@{n}"), t)
}

pub fn terminate((n, mut t): (String, T)) -> (String, T) {
    t.push(0);
    (format!("{n}$"), t)
}

pub fn variants(t: (String, T)) -> Vec<(String, T)> {
    vec![
        terminate(t.clone()),
        // t.clone(),
        terminate(flip(t.clone())),
        // flip(t.clone()),
        terminate(rev(t.clone())),
        // rev(t.clone()),
        terminate(flip(rev(t.clone()))),
        // flip(rev(t.clone())),
    ]
}

pub fn u8_minimizers((n, t): (String, T), w: usize) -> (String, T) {
    let n = format!("u8mini({n}, w={w})");

    // eprintln!("t: {t:?}");

    // convert from 1234 to ACTG
    // skip the final 0 (=$), and append a 0 at the end of the minimizers instead
    let t = AsciiSeqVec::from_ascii(
        &t[..t.len() - 1]
            .iter()
            .map(|x| b"ACGT"[*x as usize - 1])
            .collect::<Vec<_>>(),
    );
    // eprintln!("t: {t:?}");

    let k = 4;
    let mut poss = vec![];
    simd_minimizers::minimizer_positions(t.as_slice(), k, w, &mut poss);

    let mut minimizers: Vec<u8> = poss
        .into_iter()
        .map(|i| t.slice(i as usize..i as usize + k).to_word() as u8)
        .collect();
    minimizers.push(0);
    (n, minimizers)
}
