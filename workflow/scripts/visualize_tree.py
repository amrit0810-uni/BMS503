#!/usr/bin/env python3
"""
Tree Visualization — pure-Python SVG, no matplotlib or PIL required.
Reads a Newick file and writes a rectangular cladogram SVG.
Uses only Bio.Phylo (already in biopython) for tree parsing.
"""

from Bio import Phylo


def _layout(tree):
    """
    Return {id(clade): [x, y]} where:
      x = cumulative branch length from root (substitutions/site)
      y = leaf rank (integer, depth-first order); internal nodes get mean of children
    """
    coords = {}

    def _set_x(clade, parent_x=0.0):
        x = parent_x + (clade.branch_length or 0.0)
        coords[id(clade)] = [x, 0.0]
        for child in clade.clades:
            _set_x(child, x)

    rank = [0]

    def _set_y(clade):
        if clade.is_terminal():
            coords[id(clade)][1] = float(rank[0])
            rank[0] += 1
        else:
            for child in clade.clades:
                _set_y(child)
            child_ys = [coords[id(c)][1] for c in clade.clades]
            coords[id(clade)][1] = sum(child_ys) / len(child_ys)

    _set_x(tree.root)
    _set_y(tree.root)
    return coords


def tree_to_svg(tree_file, output_file):
    tree   = Phylo.read(tree_file, "newick")
    coords = _layout(tree)

    leaves   = tree.get_terminals()
    n_leaves = max(len(leaves), 2)
    all_x    = [v[0] for v in coords.values()]
    max_x    = max(all_x) if max(all_x) > 0 else 1.0

    # Layout constants
    PAD_L  = 30
    PAD_R  = 230
    PAD_T  = 30
    PAD_B  = 55
    TREE_W = 560
    ROW_H  = max(16, min(26, 560 // n_leaves))
    W      = PAD_L + TREE_W + PAD_R
    H      = PAD_T + (n_leaves - 1) * ROW_H + PAD_B

    def px(x): return PAD_L + (x / max_x) * TREE_W
    def py(y): return PAD_T + y * ROW_H

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}">',
        f'  <rect width="{W}" height="{H}" fill="white"/>',
    ]

    def _draw(clade):
        cx = px(coords[id(clade)][0])
        cy = py(coords[id(clade)][1])

        # Leaf label
        if clade.is_terminal():
            label = (clade.name or "").replace("&", "&amp;").replace("<", "&lt;")
            parts.append(
                f'  <text x="{cx + 5:.1f}" y="{cy + 4:.1f}" '
                f'font-family="monospace" font-size="11" fill="#222">{label}</text>'
            )
        elif clade.confidence is not None:
            # VeryFastTree SH-like support: stored as 0–1
            c = clade.confidence
            conf = c if c <= 1.0 else c / 100.0
            if conf > 0:
                parts.append(
                    f'  <text x="{cx + 3:.1f}" y="{cy - 3:.1f}" '
                    f'font-family="sans-serif" font-size="9" fill="#888">{conf:.2f}</text>'
                )

        for child in clade.clades:
            chx = px(coords[id(child)][0])
            chy = py(coords[id(child)][1])
            # Horizontal branch: from parent x to child x, at child y
            parts.append(
                f'  <line x1="{cx:.1f}" y1="{chy:.1f}" '
                f'x2="{chx:.1f}" y2="{chy:.1f}" stroke="#444" stroke-width="1.4"/>'
            )
            _draw(child)

        # Vertical connector spanning all children
        if clade.clades:
            ys = [py(coords[id(c)][1]) for c in clade.clades]
            parts.append(
                f'  <line x1="{cx:.1f}" y1="{min(ys):.1f}" '
                f'x2="{cx:.1f}" y2="{max(ys):.1f}" stroke="#444" stroke-width="1.4"/>'
            )

    # Root stem
    rx = px(coords[id(tree.root)][0])
    ry = py(coords[id(tree.root)][1])
    parts.append(
        f'  <line x1="{PAD_L:.1f}" y1="{ry:.1f}" '
        f'x2="{rx:.1f}" y2="{ry:.1f}" stroke="#444" stroke-width="1.4"/>'
    )

    _draw(tree.root)

    # Scale bar
    scale_val = max_x / 5
    scale_px  = TREE_W / 5
    bar_y     = H - PAD_B + 22
    parts.extend([
        f'  <line x1="{PAD_L:.1f}" y1="{bar_y}" '
        f'x2="{PAD_L + scale_px:.1f}" y2="{bar_y}" stroke="#666" stroke-width="2"/>',
        f'  <line x1="{PAD_L:.1f}" y1="{bar_y - 4}" '
        f'x2="{PAD_L:.1f}" y2="{bar_y + 4}" stroke="#666" stroke-width="1.5"/>',
        f'  <line x1="{PAD_L + scale_px:.1f}" y1="{bar_y - 4}" '
        f'x2="{PAD_L + scale_px:.1f}" y2="{bar_y + 4}" stroke="#666" stroke-width="1.5"/>',
        f'  <text x="{PAD_L:.1f}" y="{bar_y + 16}" '
        f'font-family="sans-serif" font-size="10" fill="#666">'
        f'{scale_val:.5f} substitutions/site</text>',
    ])

    parts.append('</svg>')

    with open(output_file, 'w') as f:
        f.write('\n'.join(parts))

    print(f"Tree SVG written: {output_file}  ({n_leaves} leaves)")


if __name__ == "__main__":
    input_tree = snakemake.input[0]
    output_svg = snakemake.output[0]
    tree_to_svg(input_tree, output_svg)
