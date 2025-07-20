use std::collections::hash_map::Entry;

use super::*;

pub struct Test<'t> {
    t: &'t [u8],
    map: HashMap<(usize, u8), usize>,
}

impl<'t> Test<'t> {
    pub fn new(t: &'t [u8]) -> Self {
        let n = t.len();

        let mut map = HashMap::new();

        // Simply try querying all substrings and build the map with whatever we need.
        for i in 1..n {
            let p = &t[i..n - 1];
            eprintln!(
                "Inserting pattern t[{i}..] = {}",
                str::from_utf8(p).unwrap()
            );

            let mut pos = 0;
            let mut idx = 0;
            while idx < p.len() {
                if t[pos] == p[idx] {
                    eprintln!("{idx} {pos} Match {}", p[idx] as char);
                    idx += 1;
                    pos += 1;
                    continue;
                }

                match map.entry((pos, p[idx])) {
                    Entry::Occupied(e) => {
                        pos = *e.get();
                        eprintln!("{idx} {pos} Jump to pos {pos}");
                        idx += 1;
                        pos += 1;
                        assert_eq!(p[..idx], t[pos - idx..pos]);
                    }
                    Entry::Vacant(e) => {
                        // Find an occurrence of p[..=idx] that has the longest common suffix with t[..pos]+p[idx]
                        let mut best = (0, Reverse(usize::MAX));
                        for i in 0..n - idx {
                            if &t[i..=i + idx] == &p[..=idx] {
                                best = best.max((lcs(&t[..i + idx], &t[..pos]), Reverse(i + idx)));
                            }
                        }
                        pos = best.1 .0;
                        eprintln!("{idx} {pos} Jump to pos {pos} [newly inserted]");
                        // assert!(best.0 > 0);
                        e.insert(pos);
                        pos += 1;
                        idx += 1;
                    }
                }
            }
        }

        eprintln!("\nSIZE: {}\n", map.len());

        let mut l = map.iter().collect_vec();
        l.sort();

        for ((p, c), p2) in l {
            eprintln!("pos {p} + {} => {p2}", *c as char);
        }

        Self { t, map }
    }

    pub fn locate_one(&self, p: &[u8]) -> Option<usize> {
        let mut pos = 0;
        let mut idx = 0;
        while idx < p.len() {
            if self.t[pos] == p[idx] {
                idx += 1;
                pos += 1;
                continue;
            }

            pos = *self.map.get(&(pos, p[idx]))?;
        }
        Some(pos)
    }
}

fn lcs(mut a: &[u8], mut b: &[u8]) -> usize {
    let mut l = 0;

    while a.len() > 0 && a.split_off_last() == b.split_off_last() {
        l += 1;
    }

    l
}
