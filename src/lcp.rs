use sux::{
    bits::BitVec,
    traits::{AddNumBits, Select},
};

use crate::SA;

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

impl CompactLcp {
    pub fn new(plcp: Vec<u32>) -> Self {
        // Encode plcp[i] - plcp[i-1] + 1 in (reverse) unary.
        let mut bits = sux::bits::BitVec::new(0);
        let mut prev = 0;
        for &lcp in &plcp {
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
    pub fn get(&self, sa: &SA, i: usize) -> u32 {
        let idx = if i + 1 == sa.len() { sa[0] } else { sa[i + 1] };
        (self.bits.select(idx as usize).unwrap() - 2 * idx as usize - 1) as u32
    }
    /// In bytes
    pub fn space(&self) -> usize {
        mem_dbg::MemSize::mem_size(&self.bits, mem_dbg::SizeFlags::default())
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

        let clcp = CompactLcp::new(plcp);
        for i in 0..sa.len() {
            assert_eq!(lcp[i], clcp.get(&sa, i), "i {i}");
        }
    }
}
