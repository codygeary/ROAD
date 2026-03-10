#!/usr/bin/env python3
"""
trace_pattern.py
----------------
Generates target.txt from a 2D RNA blueprint (pattern.txt).

Original Perl: Ebbe S. Andersen & Cody Geary (2013–2017)
Python port: translated from trace_pattern.pl

Usage:
    python trace_pattern.py pattern.txt [sequence.txt] > target.txt

Output lines (to stdout):
    1. Design name
    2. Dot-bracket secondary structure (with pseudoknot notation)
    3. Nucleotide constraint sequence  (T → U, X → N)

If a sequence file is given as the second argument, that sequence is
mapped onto the blueprint instead of reading constraints from the grid.
"""

import sys
from trace_common import (
    read_file, parse_grid, find_5prime, count_nts, trace_strand,
    find_base_pairs, build_structure_string, NT_MATCH, KL_CHARS,
)


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: trace_pattern.py pattern.txt [sequence.txt]")

    pattern_file  = sys.argv[1]
    sequence_file = sys.argv[2] if len(sys.argv) > 2 else None

    # ── Load optional pre-existing sequence ──────────────────────────────────
    pri: list = []
    if sequence_file:
        seq_lines = read_file(sequence_file)
        for line in seq_lines:
            print(line)   # echo the sequence file to stdout (matches Perl behaviour)
            pri = list(line)

    # ── Parse blueprint ───────────────────────────────────────────────────────
    lines = read_file(pattern_file)
    m, cross, name, kl_pattern, widest, tallest = parse_grid(lines)
    kl_pat = list(kl_pattern)

    # ── Find 5′ end and count nucleotides ─────────────────────────────────────
    r0, c0, d0 = find_5prime(m)
    nt         = count_nts(m)

    # ── Trace the strand (5′ → 3′) ────────────────────────────────────────────
    inject = pri if sequence_file else None
    seq, n, _r, _c, _d, _ok = trace_strand(m, r0, c0, d0, nt, inject_seq=inject)

    # ── Find base pairs ───────────────────────────────────────────────────────
    a_arr, b_arr, p_arr = find_base_pairs(m, n, r0, c0, d0, nt, kl_pat)

    # ── Output: name ──────────────────────────────────────────────────────────
    print(name)

    # ── Output: dot-bracket structure ─────────────────────────────────────────
    dot_bracket, _structure_map = build_structure_string(a_arr, b_arr, p_arr)
    print(dot_bracket)

    # ── Output: constraint / sequence string ──────────────────────────────────
    if not sequence_file:
        # Constraint mode: emit the nucleotide letters from the blueprint
        out_chars = []
        for ch in seq:
            if NT_MATCH.match(ch):
                out_chars.append(ch)
            elif ch == 'T':
                out_chars.append('U')
            elif ch == 'X':
                out_chars.append('N')
        print(''.join(out_chars))
    else:
        # Sequence mode: emit the injected sequence
        out_chars = []
        for ch in pri:
            if NT_MATCH.match(ch):
                out_chars.append(ch)
            elif ch == 'T':
                out_chars.append('U')
        print(''.join(out_chars))


if __name__ == '__main__':
    main()
