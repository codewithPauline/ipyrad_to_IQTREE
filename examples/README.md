# Synthetic demonstration

`demo.phy` contains eight fictional DNA sequences, each 1,200 bases long.
Regenerate it with `python examples/generate_demo.py` (standard library only).
The generator mutates a shared ancestral sequence into four pairs with seed 2026.
It is a simple demonstration, not a calibrated evolutionary simulator.

`illustrative_tree.nwk` and `figures/illustrative-tree.svg` show the intended
four-pair grouping. Branch lengths are illustrative and no support values are
assigned. **This is not an IQ-TREE result or an Ambystoma research phylogeny.**
Run the root README's demo command to infer a tree from the synthetic alignment;
the inferred topology and branch lengths may differ from the illustration.

The demo exercises alignment validation and tree inference only. It does not
exercise read preprocessing, ipyrad assembly, or the original research analysis.
