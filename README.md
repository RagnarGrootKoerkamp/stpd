# Suffix Tree Path Decompositions

This is a work-in-progress implementation of STPDs, see the preprint:
https://doi.org/10.48550/arXiv.2506.14734 and also the C++ baseline at https://github.com/regindex/STPD-index.

> Ruben Becker, Davide Cenzato, Travis Gagie, Sung-Hwan Kim, Ragnar Groot Koerkamp, Giovanni Manzini, Nicola Prezza
> Compressing Suffix Trees by Path Decompositions
> arXiv, June 2025

New ideas here:
- Using the POS- variant for incremental left-to-right construction
- Storing additional metadata for fast suffix-link computation
- MEM-finding via suffix links

Warning:
- We currently use a 2-level tiered-vec for the dynamic sparse prefix array,
  which scales badly.
