#!/usr/bin/env python3
"""
trace_analysis.py
-----------------
Restrings a designed sequence back through the 2D blueprint and produces
a full visual analysis: structure diagram, sequence overlay, GU pair
highlight, repeat-sequence map, suggested design mask, and structural
barrier (topology) analysis.

Original Perl: Ebbe S. Andersen & Cody Geary (2013–2018)
Python port: translated from trace_analysis.pl

Usage:
    python trace_analysis.py pattern.txt seq.txt

Output goes to stdout (pipe to a file or append to a design output).
"""

import sys
import random
import re
from typing import List, Optional

from trace_common import (
    read_file, parse_grid, find_5prime, count_nts, trace_strand,
    find_base_pairs, build_structure_string, build_pair_map_from_structure,
    render_grid_base,
    Grid, _get, _set,
    SYMBOL_TO_UNICODE, NT_MATCH, KL_CHARS, PAIR_MARKER, BPAIR_MARKER,
)

# ── User-configurable parameters ─────────────────────────────────────────────
COMPLEMENT_WINDOW = 10   # WC complement region length to flag (P)
DUPLICATE_WINDOW  = 9    # Duplicate region length to flag (D)
KL_DELAY          = 150  # nt delay before KLs "snap closed" in topology check
PALINDROME_WINDOW = 8    # Minimum palindrome length to flag (L)

# ── Sequence complement ───────────────────────────────────────────────────────
_COMP = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}


def _rc(seq: str) -> str:
    return ''.join(_COMP.get(b, 'N') for b in reversed(seq))


# ── Repeat / pattern counter ─────────────────────────────────────────────────

def count_repeats(trial_sol: List[str], nt: int) -> dict:
    """
    Scan trial_sol for unwanted sequence patterns.
    Returns a dict with repeat_map list and counts.

    Note: this version uses 8-nt runs (vs dragon.py which uses 6-nt runs) to
    match the original trace_analysis.pl thresholds.
    """
    repeat_map    = ['-'] * nt
    seq_str       = ''.join(trial_sol[:nt])

    # Complementary regions
    for i in range(nt - COMPLEMENT_WINDOW):
        window = seq_str[i:i + COMPLEMENT_WINDOW]
        rc     = _rc(window)
        pos    = 0
        while True:
            found = seq_str.find(rc, pos)
            if found == -1:
                break
            for k in range(COMPLEMENT_WINDOW):
                if i + k < nt:     repeat_map[i + k]     = 'P'
                if found + k < nt: repeat_map[found + k] = 'P'
            pos = found + 1

    complement_zones = sum(1 for c in repeat_map if c == 'P')

    # Palindromes
    for i in range(nt - PALINDROME_WINDOW):
        window = seq_str[i:i + PALINDROME_WINDOW]
        rc     = _rc(window)
        pos    = i
        while pos < nt:
            p1 = seq_str.find(rc, pos)
            p2 = seq_str.find(window, pos)
            if p1 == -1 or p2 == -1:
                break
            if p1 == p2:
                for k in range(PALINDROME_WINDOW):
                    if i + k < nt:  repeat_map[i + k]  = 'L'
                    if p1 + k < nt: repeat_map[p1 + k] = 'L'
            pos = max(p1, p2) + 1

    palindromes = sum(1 for c in repeat_map if c == 'L')

    # Duplicated regions
    for i in range(nt - DUPLICATE_WINDOW):
        window = seq_str[i:i + DUPLICATE_WINDOW]
        pos    = i
        while True:
            found = seq_str.find(window, pos)
            if found == -1:
                break
            if found > i:
                for k in range(DUPLICATE_WINDOW):
                    if i + k < nt:     repeat_map[i + k]     = 'D'
                    if found + k < nt: repeat_map[found + k] = 'D'
            pos = found + 1

    duplicate_zones = sum(1 for c in repeat_map if c == 'D')

    # AU-rich / GC-rich runs of 8+
    pattern_repeats = 0
    for i in range(nt - 7):
        window8 = trial_sol[i:i + 8]
        if all(b in ('A', 'U') for b in window8):
            pattern_repeats += 1
            for k in range(8):
                repeat_map[i + k] = 'W'
        if all(b in ('C', 'G') for b in window8):
            pattern_repeats += 1
            for k in range(8):
                repeat_map[i + k] = 'S'

    # Poly-homopolymer runs of 4+
    poly_repeats = 0
    for i in range(nt - 6):
        if all(trial_sol[i + k] == trial_sol[i] for k in range(1, 5)):
            poly_repeats += 1
            base = trial_sol[i]
            for k in range(5):
                if i + k < nt:
                    repeat_map[i + k] = base

    # Restriction sites
    _SITES = [
        'GGUCUC', 'GAGACC', 'GAAGAC', 'GUCUUC',
        'CGUCUC', 'GAGACG',
        'GCUCUUC', 'GAAGAGC', 'AUCUGUU',
    ]
    restriction_sites = 0
    for i in range(nt - 3):
        for site in _SITES:
            end = i + len(site)
            if end <= nt and seq_str[i:end] == site:
                restriction_sites += 1
                for k in range(len(site)):
                    repeat_map[i + k] = 'X'
                break

    return dict(
        repeat_map=repeat_map,
        complement_zones=complement_zones,
        palindromes=palindromes,
        duplicate_zones=duplicate_zones,
        pattern_repeats=pattern_repeats,
        poly_repeats=poly_repeats,
        restriction_sites=restriction_sites,
    )


# ── Grid rendering helpers ────────────────────────────────────────────────────

def _render_skeleton(m: Grid, strand_dir: Grid, cross: Grid,
                     node_fn=None) -> str:
    """
    Render the grid skeleton (connectors + optional node_fn for NT cells).
    node_fn(i, j) → character string for nucleotide cells; None → standard.
    """
    out = []
    for i in range(1000):
        row = m.get(i)
        if row is None:
            break
        for j in sorted(row.keys()):
            ch = row[j]
            if not ch:
                continue

            if ch == '\n':   out.append('\n');                     continue
            if ch == ' ':    out.append(' ');                      continue
            if ch in ('3', '5'): out.append(ch);                  continue
            if ch == 'T':    out.append('─');                      continue

            # Box-drawing
            sym = SYMBOL_TO_UNICODE.get(ch)
            if sym:          out.append(sym);                      continue

            if ch == 'b':    out.append('─');                      continue
            if ch == 'i':
                if cross and _get(cross, i, j) == '^':
                    out.append('^')
                else:
                    out.append('│')
                continue
            if ch == 'x':    out.append('┼');                      continue

            # Pair connectors — direction-dependent rendering
            if ch in ('p', '!', '*'):
                sd = _get(strand_dir, i, j)
                if sd in ('down', 'up'):
                    # strand runs vertically → connector is horizontal '─'
                    out.append('─')
                else:
                    if ch == 'p':
                        out.append('┊')
                    elif ch == '!':
                        out.append('!')
                    else:
                        out.append('*')
                continue

            # Nucleotide cells
            if NT_MATCH.match(ch):
                if node_fn:
                    out.append(node_fn(i, j))
                else:
                    out.append(ch)
                continue

            out.append(ch)
    return ''.join(out)


def _strand_symbol(m: Grid, strand_dir: Grid, i: int, j: int) -> str:
    """Return the strand-path symbol for a nucleotide cell."""
    sd = _get(strand_dir, i, j)
    if sd in ('down', 'up'):
        return '│'
    return '─'


def _sequence_symbol(m: Grid, n: Grid, seq_arr: List[str], i: int, j: int) -> str:
    pos = _get(n, i, j)
    if pos:
        return seq_arr[pos - 1]
    ch = _get(m, i, j)
    return ch if NT_MATCH.match(ch) else ch


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 3:
        sys.exit("Usage: trace_analysis.py pattern.txt seq.txt")

    pattern_file  = sys.argv[1]
    sequence_file = sys.argv[2]

    # ── Load sequence file ────────────────────────────────────────────────────
    seq_lines = read_file(sequence_file)
    pri: List[str] = []
    for line in seq_lines:
        print(line)
        pri = list(line)

    # ── Parse blueprint ───────────────────────────────────────────────────────
    lines = read_file(pattern_file)
    m, cross, name, kl_pattern, widest, tallest = parse_grid(lines)
    kl_pat = list(kl_pattern)

    # strand_dir tracks the direction of the strand at each cell (for rendering)
    strand_dir: Grid = {}

    # ── Find 5′ end ───────────────────────────────────────────────────────────
    r0, c0, d0 = find_5prime(m)
    _set(strand_dir, r0, c0, d0)
    nt = count_nts(m)

    # ── Trace the strand, injecting the sequence ──────────────────────────────
    seq_arr, n, _r, _c, _d, _ok = trace_strand(
        m, r0, c0, d0, nt, inject_seq=pri
    )
    # Populate strand_dir during a second tracing pass (we need both seq and dirs)
    _populate_strand_dir(m, strand_dir, r0, c0, d0, nt)

    # ── Find base pairs ───────────────────────────────────────────────────────
    a_arr, b_arr, p_arr = find_base_pairs(
        m, n, r0, c0, d0, nt, kl_pat, strand_dir=strand_dir
    )

    # ── Print: name, structure, sequence ─────────────────────────────────────
    print(name)

    dot_bracket, structure_map = build_structure_string(a_arr, b_arr, p_arr)
    print(dot_bracket)

    sequence = ''
    for ch in (pri if sequence_file else seq_arr):
        if NT_MATCH.match(ch):
            sequence += ch
        elif ch == 'T':
            sequence += 'U'
        elif ch == 'X':
            sequence += 'N'
    print(sequence)

    print('\n\n')
    print('My Structure map:  ')
    print(structure_map, '\n\n')

    pair_map = build_pair_map_from_structure(structure_map)

    # ── Scrub non-AUCG positions in the sequence ──────────────────────────────
    sequence_array = list(sequence)
    _scrub_sequence(sequence_array, pair_map, nt)

    # ── 2D diagram with sequence ──────────────────────────────────────────────
    print('\n2D diagram with sequence')
    _print_2d_with_sequence(m, n, strand_dir, cross, seq_arr, a_arr, b_arr)

    # ── Strand path ───────────────────────────────────────────────────────────
    print('\n\nStrand Path')
    out = _render_skeleton(m, strand_dir, cross,
                           node_fn=lambda i, j: _strand_symbol(m, strand_dir, i, j))
    print(out)

    # ── Topology / structural barriers check ──────────────────────────────────
    print('\n\nHighlighting Structural Barriers\n\n')
    struc_array = list(structure_map)
    barriers    = _compute_barriers(struc_array, pair_map, nt)

    out = _render_skeleton(m, strand_dir, cross,
        node_fn=lambda i, j: barriers[(_get(n, i, j) or 1) - 1])
    print(out)

    # ── GU pair highlighter ───────────────────────────────────────────────────
    print('\n\nHighlighting GU Pairs')
    wobble_sym = '\u25e6'   # ◦  (hollow bullet, used for non-GU)
    wobbles_seq = [wobble_sym] * nt
    for i in range(nt):
        ai = sequence_array[i]
        mi = pair_map[i]
        bi = sequence_array[mi] if mi != i else ''
        if {ai, bi} == {'G', 'U'} or (ai == 'K' and bi == 'K'):
            wobbles_seq[i]  = ai
            wobbles_seq[mi] = bi

    out = _render_skeleton(m, strand_dir, cross,
        node_fn=lambda i, j: wobbles_seq[(_get(n, i, j) or 1) - 1])
    print(out)

    # ── Repeat sequence analysis ──────────────────────────────────────────────
    sequence_check = all(b in ('A', 'U', 'C', 'G') for b in sequence_array[:nt])

    if not sequence_check:
        print('\nSequence not suitable for pattern search. '
              'Must only contain A, U, C and G.')
        return

    rpt = count_repeats(sequence_array[:nt], nt)
    repeat_map     = rpt['repeat_map']
    comp_zones     = rpt['complement_zones']
    pals           = rpt['palindromes']
    dup_zones      = rpt['duplicate_zones']
    pat_repeats    = rpt['pattern_repeats']
    poly_repeats   = rpt['poly_repeats']
    restrict_sites = rpt['restriction_sites']

    print('\n\nHighlighting Repeat Sequences \n\n')
    print(f' WC complement region (P) {COMPLEMENT_WINDOW} or longer: {comp_zones} nts')
    print(f' Palindromes {PALINDROME_WINDOW} or longer (L): {pals}')
    print(f' Duplicated region (D) {DUPLICATE_WINDOW} or longer: {dup_zones} nts')
    print(f' Strong/Weak region (S/W) 8 or longer: {pat_repeats} nts')
    print(f' 5 or more in a row of the same nucleotide (A,U,C,G): {poly_repeats} nts')
    print(f' Common restriction site (X): {restrict_sites} nts')

    # Repeat map overlay
    hollow = '\u25e6'
    def _repeat_node(i, j):
        pos = _get(n, i, j)
        if not pos:
            return _get(m, i, j)
        rm = repeat_map[pos - 1]
        return hollow if rm == '-' else rm

    out = _render_skeleton(m, strand_dir, cross, node_fn=_repeat_node)
    print(out)

    # ── Suggested design mask ─────────────────────────────────────────────────
    print(' Printing a suggested design mask')
    print(' Sites to be replaced are marked (N) as are their partners. ')

    def _mask_node(i, j):
        pos = _get(n, i, j)
        if not pos:
            ch = _get(m, i, j)
            return ch if NT_MATCH.match(ch) else ch
        rm_i  = repeat_map[pos - 1]
        rm_mp = repeat_map[pair_map[pos - 1]]
        if rm_i == '-' and rm_mp == '-':
            return _get(m, i, j)
        return 'N'

    out = _render_skeleton(m, strand_dir, cross, node_fn=_mask_node)
    print(out)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _populate_strand_dir(m: Grid, strand_dir: Grid,
                          r0: int, c0: int, d0: str, nt: int) -> None:
    """
    Walk the strand again purely to record direction at every cell.
    (A lightweight re-trace that doesn't collect sequence data.)
    """
    from trace_common import HORIZ_MOVE, VERT_MOVE, HORIZ_STEP, VERT_STEP, _apply_crossover, DIGIT

    r, c, d = r0, c0, d0
    _set(strand_dir, r, c, d)

    for _ in range(nt + 10000):
        _set(strand_dir, r, c, d)

        moved = False
        if d == 'right' and HORIZ_MOVE.match(_get(m, r, c + 1)):
            c += 1; moved = True
        elif d == 'left' and HORIZ_MOVE.match(_get(m, r, c - 1)):
            c -= 1; moved = True
        elif d == 'up' and VERT_MOVE.match(_get(m, r - 1, c)):
            r -= 1; moved = True
        elif d == 'down' and VERT_MOVE.match(_get(m, r + 1, c)):
            r += 1; moved = True
        else:
            if d == 'right' and HORIZ_STEP.match(_get(m, r, c + 1)):
                c += 1
            elif d == 'left' and HORIZ_STEP.match(_get(m, r, c - 1)):
                c -= 1
            elif d == 'down' and VERT_STEP.match(_get(m, r + 1, c)):
                r += 1
            elif d == 'up' and VERT_STEP.match(_get(m, r - 1, c)):
                r -= 1

        d, r, c = _apply_crossover(m, d, r, c)
        _set(strand_dir, r, c, d)

        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nb = _get(m, r + dr, c + dc)
            if DIGIT.match(nb) and int(nb) == 3:
                return


def _scrub_sequence(sequence_array: List[str], pair_map: List[int], nt: int) -> None:
    """Replace non-AUCG positions with random Watson-Crick pairs in-place."""
    bases = ['A', 'U', 'G', 'C']
    pairs_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    for i in range(nt):
        if sequence_array[i] not in ('A', 'U', 'C', 'G'):
            pick = random.randint(0, 3)
            b    = bases[pick]
            sequence_array[i]           = b
            sequence_array[pair_map[i]] = pairs_map[b]


def _print_2d_with_sequence(m: Grid, n: Grid, strand_dir: Grid, cross: Grid,
                             seq_arr: List[str], a_arr, b_arr) -> None:
    """Render the 2D diagram with the actual sequence and WC pair checking."""
    out = []
    for i in range(1000):
        row = m.get(i)
        if row is None:
            break
        for j in sorted(row.keys()):
            ch = row.get(j, '')
            if not ch:
                continue
            if ch == '\n':     out.append('\n');          continue
            if ch == ' ':      out.append(' ');           continue
            if ch in ('3','5'):out.append(ch);            continue
            if ch == 'T':      out.append('U');           continue

            from trace_common import SYMBOL_TO_UNICODE
            sym = SYMBOL_TO_UNICODE.get(ch)
            if sym:            out.append(sym);           continue

            if ch == 'b':
                # Horizontal pair connector — check WC validity
                left  = _get(m, i, j - 1)
                right = _get(m, i, j + 1)
                if _is_valid_pair(left, right):
                    out.append('=')
                else:
                    out.append('?')
                continue

            if ch == 'i':
                from trace_common import _get as cget
                if _get(cross, i, j) == '^':
                    out.append('^')
                else:
                    out.append('│')
                continue

            if ch == 'x':      out.append('┼');           continue

            if ch == 'p':
                # Vertical pair connector — check WC validity
                above = _get(m, i - 1, j)
                below = _get(m, i + 1, j)
                if _is_valid_pair(above, below):
                    out.append('┊')
                else:
                    out.append('?')
                continue

            if ch == '!':      out.append('!');           continue

            if NT_MATCH.match(ch):
                out.append(ch)
                continue

            out.append(ch)

    print(''.join(out))


def _is_valid_pair(a: str, b: str) -> bool:
    valid = {frozenset({'G', 'C'}), frozenset({'A', 'U'}), frozenset({'G', 'U'})}
    ambig = {'N', 'S', 'K', 'R', 'Y', 'M', 'W', 'V', 'H', 'B', 'D'}
    if a in ambig or b in ambig:
        return True
    return frozenset({a, b}) in valid


def _compute_barriers(struc_array: List[str], pair_map: List[int], nt: int) -> List[str]:
    """
    Topology checker: marks structural barriers in the folding pathway.

    x = open pair (potential barrier)
    X = blocked pair (barrier reached)
    ~ = ambiguous / partially blocked
    ◦ = clear
    """
    hollow     = '\u25e6'
    barriers   = [hollow] * nt
    topo_count = 0

    for i in range(nt):
        ch = struc_array[i]
        mi = pair_map[i]

        if ch == '.':
            barriers[i] = hollow
            topo_count  = 0

        elif ch == '(':
            barriers[i] = 'x'
            topo_count  = 0

        elif ch == ')':
            if barriers[mi] == 'x':
                barriers[i]  = hollow
                barriers[mi] = hollow
                topo_count   = 0
            elif barriers[mi] == 'X':
                if topo_count > 5:
                    barriers[i]  = 'X'
                    barriers[mi] = '~'
                else:
                    barriers[i]  = '~'
                    barriers[mi] = '~'
                topo_count += 1

        # KL closure: after KL_DELAY nts, mark intervening sequence as blocked
        if i > KL_DELAY:
            delayed_i = i - KL_DELAY
            if struc_array[delayed_i] == ']':
                barriers[delayed_i]         = hollow
                barriers[pair_map[delayed_i]] = hollow
                for k in range(pair_map[delayed_i], delayed_i):
                    if barriers[k] == 'x':
                        barriers[k] = 'X'

    return barriers


if __name__ == '__main__':
    main()
