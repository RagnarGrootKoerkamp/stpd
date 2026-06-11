use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};
use sux::{
    bits::BitVec,
    dict::{EfSeq, EliasFanoBuilder},
    traits::{AddNumBits, IndexedSeq, Select},
};

use crate::SA;

pub trait Lcp: AsRef<Self> + Sync + serde::Serialize + for<'a> serde::Deserialize<'a> {
    fn new(sa: &SA, plcp: &Vec<u32>) -> Self;
    fn get(&self, sa: &SA, i: usize) -> u32;
    fn space(&self) -> usize;
}

/// See https://ae.iti.kit.edu/download/ti_lec11_1.pdf
#[derive(serde::Serialize, serde::Deserialize)]
pub struct PlainLcp {
    lcp: Vec<u32>,
}

impl AsRef<PlainLcp> for PlainLcp {
    fn as_ref(&self) -> &PlainLcp {
        self
    }
}

impl Lcp for PlainLcp {
    fn new(sa: &SA, plcp: &Vec<u32>) -> Self {
        Self {
            lcp: (0..sa.len())
                .into_par_iter()
                .map(|i| plcp[sa[if i + 1 < sa.len() { i + 1 } else { 0 }] as usize])
                .collect(),
        }
    }
    fn get(&self, _sa: &SA, i: usize) -> u32 {
        self.lcp[i]
    }
    /// In bytes
    fn space(&self) -> usize {
        mem_dbg::MemSize::mem_size(&self.lcp, mem_dbg::SizeFlags::default())
    }
}

/// See https://ae.iti.kit.edu/download/ti_lec11_1.pdf
#[derive(serde::Serialize, serde::Deserialize)]
pub struct CompactLcp {
    bits: sux::rank_sel::SelectAdaptConst<AddNumBits<BitVec>>,
}

impl AsRef<CompactLcp> for CompactLcp {
    fn as_ref(&self) -> &CompactLcp {
        self
    }
}

impl Lcp for CompactLcp {
    fn new(_sa: &SA, plcp: &Vec<u32>) -> Self {
        // Encode plcp[i] - plcp[i-1] + 1 in (reverse) unary.
        let mut bits = sux::bits::BitVec::new(0);
        let mut prev = 0;
        for &lcp in plcp {
            let delta = (lcp + 1).strict_sub(prev);
            bits.extend(std::iter::repeat(false).take(delta as usize));
            bits.push(true);
            prev = lcp;
        }
        bits.push(true);
        Self {
            bits: sux::rank_sel::SelectAdaptConst::new(bits.into()),
        }
    }
    fn get(&self, sa: &SA, i: usize) -> u32 {
        let idx = if i + 1 == sa.len() { sa[0] } else { sa[i + 1] };
        (self.bits.select(idx as usize).unwrap() - 2 * idx as usize - 1) as u32
    }
    /// In bytes
    fn space(&self) -> usize {
        mem_dbg::MemSize::mem_size(&self.bits, mem_dbg::SizeFlags::default())
    }
}

/// See https://ae.iti.kit.edu/download/ti_lec11_1.pdf
#[derive(serde::Serialize, serde::Deserialize)]
pub struct EfLcp {
    ef: EfSeq<u64>,
}

impl AsRef<EfLcp> for EfLcp {
    fn as_ref(&self) -> &EfLcp {
        self
    }
}

impl Lcp for EfLcp {
    fn new(sa: &SA, plcp: &Vec<u32>) -> Self {
        let n = sa.len();
        let sum = plcp.iter().map(|x| *x as u64).sum();
        eprintln!("Sum of PLCP is {sum}");
        eprintln!("avg PLCP is {}", sum as f64 / n as f64);
        let mut builder = EliasFanoBuilder::new(n + 1, sum);
        let mut last = 0;
        builder.push(last);
        eprintln!("Pushing EfLcp ..");
        let mut buf = vec![];
        // Split data into chunks of size 2^32 that we collect in parallel.
        for start in (0..n).step_by(1 << 25) {
            let end = (start + (1 << 25)).min(n);

            (start..end)
                .into_par_iter()
                .map(|i| {
                    if i + 32 < n {
                        // prefetch ahead.
                        let i = i + 32;
                        let idx = if i + 1 == sa.len() { sa[0] } else { sa[i + 1] };
                        prefetch_index::prefetch_index(plcp.as_slice(), idx as usize);
                    }
                    let idx = if i + 1 == sa.len() { sa[0] } else { sa[i + 1] };
                    plcp[idx as usize] as u64
                })
                .collect_into_vec(&mut buf);
            for &lcp in &buf {
                last += lcp;
                builder.push(last);
            }
            buf.clear();
        }
        drop(buf);
        eprintln!("Building EfLcp ..");
        Self {
            ef: builder.build_with_seq(),
        }
    }
    fn get(&self, _sa: &SA, i: usize) -> u32 {
        (self.ef.get(i + 1) - self.ef.get(i)) as u32
    }
    /// In bytes
    fn space(&self) -> usize {
        mem_dbg::MemSize::mem_size(&self.ef, mem_dbg::SizeFlags::default())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::LCP;

    #[test]
    fn test_compact_lcp() {
        let t = b"bananaannabannabancbabanabacncbcncnncnabana";

        let sa: SA;
        let plcp: LCP;
        let mut lcp: LCP;
        {
            let sa_builder = libsais::SuffixArrayConstruction::for_text(t.as_slice())
                .in_owned_buffer64()
                .multi_threaded(libsais::ThreadCount::openmp_default())
                .run()
                .unwrap();

            let plcp_builder = sa_builder
                .plcp_construction()
                .multi_threaded(libsais::ThreadCount::openmp_default())
                .run()
                .unwrap();

            let (_sa, _plcp, _) = plcp_builder.into_parts();
            sa = _sa.into_iter().map(|x| x as u32).collect();
            plcp = _plcp.into_iter().map(|x| x as u32).collect();

            let n = t.len();
            lcp = (0..n).map(|i| plcp[sa[i] as usize]).collect();

            lcp.remove(0);
            lcp.push(0);
        }
        eprintln!("sa   {sa:?}");
        eprintln!("plcp {plcp:?}");
        eprintln!("lcp  {lcp:?}");

        let clcp = CompactLcp::new(&sa, &plcp);
        for i in 0..sa.len() {
            assert_eq!(lcp[i], clcp.get(&sa, i), "i {i}");
        }
    }
}
