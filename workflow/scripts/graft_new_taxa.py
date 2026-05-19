#!/usr/bin/env python3
"""
Graft taxa that appear in an alignment but not in an existing Newick tree onto
the root of that tree.  The resulting tree is written to disk for use as a
FastTree --intree starting topology; FastTree then optimises branch lengths and
placement for all sequences, including the newly grafted ones.

Usage:
    python3 graft_new_taxa.py <existing_tree.nwk> <alignment.fasta> <output_tree.nwk>
"""
import sys
import shutil
from Bio import Phylo, SeqIO
from Bio.Phylo.BaseTree import Clade

existing_tree_file = sys.argv[1]
alignment_file = sys.argv[2]
output_tree_file = sys.argv[3]

trees = list(Phylo.parse(existing_tree_file, "newick"))
tree = trees[0]

existing_taxa = {c.name for c in tree.get_terminals() if c.name}
all_taxa = [r.id for r in SeqIO.parse(alignment_file, "fasta")]
all_taxa_set = set(all_taxa)
new_taxa = [t for t in all_taxa if t not in existing_taxa]
stale_taxa = [t for t in existing_taxa if t not in all_taxa_set]

if stale_taxa:
    print(f"Pruning {len(stale_taxa)} stale taxa from master tree: {stale_taxa}")
    for taxon in stale_taxa:
        tree.prune(taxon)

if new_taxa:
    # FastTree --intree requires a strictly binary tree; adding new taxa as root
    # children creates a polytomy that causes an assertion failure.  Signal the
    # caller to build from scratch instead.
    print(f"New taxa detected ({len(new_taxa)}): {new_taxa} — caller should build tree from scratch")
    if stale_taxa:
        # Write the pruned tree so the caller still has a valid starting point
        Phylo.write(tree, output_tree_file, "newick")
    sys.exit(1)

if stale_taxa:
    Phylo.write(tree, output_tree_file, "newick")
else:
    print(f"No new or stale taxa — copying master tree as starting topology ({len(existing_taxa)} taxa)")
    shutil.copy(existing_tree_file, output_tree_file)
