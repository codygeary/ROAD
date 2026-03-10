"""
trace_common.py
---------------
Shared grid-parsing, strand-tracing, and base-pair-finding logic used by
both trace_pattern.py and trace_analysis.py.

Original Perl: Ebbe S. Andersen & Cody Geary (2013–2018)
Python port: translated from trace_pattern.pl / trace_analysis.pl

Blueprint notation:
  5          → 5′ end marker
  3          → 3′ end marker
  N/X        → unconstrained nucleotide position
  A/U/C/G    → constrained nucleotide
  R/Y/K/M/S/W/V/H/B/D → IUPAC ambiguity codes
  T          → thymine (converted to U)
  -  or ─   → horizontal strand segment
  |  or │   → vertical strand segment (i)
  :  or ┊   → base-pair connector (p)
  +  or ┼   → crossover (x)
  !          → extra/forced base pair
  *          → kissing-loop base pair
  =  or b   → horizontal base-pair connector (b)
  ^          → annotated vertical
  #          → comment: skip rest of line
  >name      → design name
  @pattern   → KL matching pattern

ASCII art corner shortcuts (context-sensitive):
  \\ + NT-char  → L  (curve: strand turns right going up / left going down)
  \\ + other   → 7  (curve: strand turns left going up)
  /  + NT-char  → r  (curve: strand turns right going down)
  /  + other   → J  (curve: strand turns left going up from below)

Unicode box-drawing equivalents are also accepted.
"""

import sys
import re
from typing import Dict, List, Optional, Tuple

# ── Unicode box-drawing → internal symbol ────────────────────────────────────
# The Perl matched individual bytes of UTF-8 sequences; here we match the
# proper Python Unicode characters instead.

_UNICODE_TO_SYMBOL: Dict[str, str] = {
    '─': '-',   # U+2500  BOX DRAWINGS LIGHT HORIZONTAL
    '│': 'i',   # U+2502  BOX DRAWINGS LIGHT VERTICAL
    '┊': 'p',   # U+250A  BOX DRAWINGS LIGHT DOUBLE DASH VERTICAL
    '┼': 'x',   # U+253C  BOX DRAWINGS LIGHT VERTICAL AND HORIZONTAL
    '╭': 'r',   # U+256D  BOX DRAWINGS LIGHT ARC DOWN AND RIGHT
    '╮': '7',   # U+256E  BOX DRAWINGS LIGHT ARC DOWN AND LEFT
    '╯': 'J',   # U+256F  BOX DRAWINGS LIGHT ARC UP AND LEFT
    '╰': 'L',   # U+2570  BOX DRAWINGS LIGHT ARC UP AND RIGHT
    '\u00a0': ' ',  # non-breaking space → plain space
}

# Internal symbol → Unicode for output rendering
SYMBOL_TO_UNICODE: Dict[str, str] = {
    '-': '─', 'i': '│', 'p': '┊', 'x': '┼',
    'r': '╭', '7': '╮', 'J': '╯', 'L': '╰',
}

# Characters that count as nucleotide positions on the blueprint
NT_CHARS = re.compile(r'[NXGACURYKMSWVHBDT]')

# Characters that the strand tracer can move through (horizontal)
HORIZ_MOVE   = re.compile(r'[xNXGACURYKMSWVHBDT\-]')
VERT_MOVE    = re.compile(r'[xNXGACURYKMSWVHBDTi]')
NT_MATCH     = re.compile(r'[NXGACURYKMSWVHBDT]')
HORIZ_STEP   = re.compile(r'[x\-]')
VERT_STEP    = re.compile(r'[xi]')
PAIR_MARKER  = re.compile(r'[!p*]')
BPAIR_MARKER = re.compile(r'[b*]')
DIGIT        = re.compile(r'\d')
KL_CHARS     = re.compile(r'[ABCDEFGHIJ1234567890]')


# ── File reading ──────────────────────────────────────────────────────────────

def read_file(filename: str) -> List[str]:
    """Read a file, normalise line endings, return list of lines."""
    with open(filename, 'r', encoding='utf-8', errors='replace') as f:
        content = f.read()
    content = content.replace('\r\n', '\n').replace('\r', '\n')
    return content.split('\n')


# ── Grid types ────────────────────────────────────────────────────────────────

# The grid is a dict-of-dicts: m[row][col] = symbol_char
Grid = Dict[int, Dict[int, str]]


def _get(m: Grid, r: int, c: int, default: str = '') -> str:
    return m.get(r, {}).get(c, default)


def _set(m: Grid, r: int, c: int, val: str) -> None:
    if r not in m:
        m[r] = {}
    m[r][c] = val


# ── Grid parser ───────────────────────────────────────────────────────────────

def parse_grid(lines: List[str]) -> Tuple[Grid, Grid, str, str, int, int]:
    """
    Parse the blueprint lines into the internal 2D grid.

    Returns:
        m          – the symbol grid
        cross      – secondary grid for '^' annotations
        name       – design name extracted from '>' line
        kl_pattern – KL pattern string from '@' line
        widest     – rightmost active column
        tallest    – number of rows
    """
    m: Grid     = {}
    cross: Grid = {}
    name        = 'Untitled'
    kl_pattern  = ''
    maxlength   = 0
    widest      = 0
    j           = 0   # row index

    for line in lines:
        l     = len(line)
        k     = 0       # column index within this row
        wide  = 0       # rightmost active column in this row

        if j == 0:
            maxlength = l

        i = 0
        while i < l:
            cols = line[i]

            # Translate Unicode box-drawing characters first
            cols = _UNICODE_TO_SYMBOL.get(cols, cols)

            # Comments and metadata
            if cols == '\t':
                sys.exit("Error: Pattern file contains tabs. Replace with spaces.\n")
            if cols == '#':
                break   # rest of line is a comment
            if cols == '>':
                name = line[i + 1:i + 33].strip()
                break
            if cols == '@':
                kl_pattern = line[i + 1:i + 33].strip()
                break

            # Add left-edge space on first character of each row
            if k == 0:
                _set(m, j, k, ' ')
                k += 1

            # ── Map input characters to internal symbols ──────────────
            next_ch = line[i + 1] if i + 1 < l else ''
            nt_next = bool(re.match(r'[+xNXGACURYKMSWVHBDT35\-]', next_ch))

            placed = True
            if   cols == ' ':             sym = ' '
            elif cols == '*':             sym = '*';  wide = k
            elif cols == '/':
                sym = 'r' if nt_next else 'J';        wide = k
            elif cols == '\\':
                sym = 'L' if nt_next else '7';        wide = k
            elif cols in ('-', '─'):      sym = '-';  wide = k
            elif cols in ('|', '│'):      sym = 'i';  wide = k
            elif cols == '^':
                sym = 'i'
                _set(cross, j, k, '^')
                wide = k
            elif cols in (':', ';', '┊'): sym = 'p';  wide = k
            elif cols == '!':             sym = '!';  wide = k
            elif cols in ('+', '┼'):      sym = 'x';  wide = k
            elif cols == '=':             sym = 'b';  wide = k
            # Already-translated symbols from unicode map
            elif cols in ('r', '7', 'J', 'L'):
                sym = cols; wide = k
            elif re.match(r'\w', cols):   sym = cols; wide = k
            else:
                placed = False

            if placed:
                _set(m, j, k, sym)
                k += 1

            # If this row is now wider than any previous, pad all previous rows
            if k > maxlength:
                for mi in range(j + 1):
                    for pn in range(maxlength, k):
                        existing = _get(m, mi, pn)
                        if existing == '\n':
                            _set(m, mi, pn, ' ')
                            _set(m, mi, pn + 1, '\n')
                        elif existing == '':
                            _set(m, mi, pn, 'U')
                            _set(m, mi, pn + 1, '\n')
                maxlength = k

            if wide > widest:
                widest = wide

            i += 1

        # Pad short rows to maxlength
        endpadding = maxlength - k
        for _ in range(endpadding):
            _set(m, j, k, ' ')
            k += 1
        _set(m, j, k, '\n')

        j += 1

    tallest = j

    # Trim excess whitespace on the right edge
    for ri in range(l if lines else 0):
        for ci in range(widest + 1, maxlength):
            _set(m, ri, ci, '')
        _set(m, widest, ci, '\n')

    # Add two blank rows at the bottom
    for ri in [tallest, tallest + 1]:
        for ci in range(widest + 2):
            _set(m, ri, ci, ' ')
        _set(m, ri, widest + 2, '\n')

    return m, cross, name, kl_pattern, widest, tallest


# ── 5′ end finder ─────────────────────────────────────────────────────────────

def find_5prime(m: Grid) -> Tuple[int, int, str]:
    """
    Scan the grid for the cell labelled '5' and determine the initial
    tracing direction from its neighbours.

    Returns (row, col, direction) where direction ∈ {'right','left','up','down'}.
    """
    for i in range(1000):
        for j in range(1000):
            cell = _get(m, i, j)
            if cell and DIGIT.match(cell) and int(cell) == 5:
                d = 'left'   # default
                if re.match(r'[NXGACURYKMSWVHBDT\-]', _get(m, i, j + 1)):
                    d = 'right'
                if re.match(r'[NXGACURYKMSWVHBDT\-]', _get(m, i, j - 1)):
                    d = 'left'
                if re.match(r'[NXGACURYKMSWVHBDTi]',  _get(m, i + 1, j)):
                    d = 'down'
                if re.match(r'[NXGACURYKMSWVHBDTi]',  _get(m, i - 1, j)):
                    d = 'up'
                return i, j, d
    raise RuntimeError("5′ end not found in blueprint. Is '5' present?")


# ── Nucleotide counter ────────────────────────────────────────────────────────

def count_nts(m: Grid) -> int:
    """Count all cells containing a nucleotide character."""
    total = 0
    for i in range(1000):
        for j in range(1000):
            if NT_MATCH.match(_get(m, i, j)):
                total += 1
    return total


# ── Strand tracer ─────────────────────────────────────────────────────────────

def trace_strand(
    m: Grid, r0: int, c0: int, d0: str, nt: int,
    inject_seq: Optional[List[str]] = None
) -> Tuple[List[str], Grid, int, int, str, bool]:
    """
    Walk the strand from 5′ to 3′, collecting the nucleotide sequence and
    recording sequential position numbers in the n-grid.

    If inject_seq is provided (trace_analysis mode), nucleotide cells are
    overwritten with the characters from inject_seq instead of being read.

    Returns:
        seq        – list of nucleotide characters in 5′→3′ order
        n          – position-number grid (n[r][c] = 1-based nt index)
        r_end, c_end, d_end  – position/direction when 3′ end was reached
        success    – True if the 3′ marker (value==3) was found
    """
    m_work = m   # we may mutate m in inject mode
    n: Grid = {}
    seq: List[str] = []
    num  = 0
    r, c = r0, c0
    d    = d0
    success = False

    for _ in range(nt + 10000):
        # ── Advance in current direction ──────────────────────────────
        if d == 'right' and HORIZ_MOVE.match(_get(m_work, r, c + 1)):
            if NT_MATCH.match(_get(m_work, r, c + 1)):
                if inject_seq is not None:
                    _set(m_work, r, c + 1, inject_seq[num])
                num += 1
                _set(n, r, c + 1, num)
                seq.append(_get(m_work, r, c + 1))
            c += 1

        elif d == 'left' and HORIZ_MOVE.match(_get(m_work, r, c - 1)):
            if NT_MATCH.match(_get(m_work, r, c - 1)):
                if inject_seq is not None:
                    _set(m_work, r, c - 1, inject_seq[num])
                num += 1
                _set(n, r, c - 1, num)
                seq.append(_get(m_work, r, c - 1))
            c -= 1

        elif d == 'up' and VERT_MOVE.match(_get(m_work, r - 1, c)):
            if NT_MATCH.match(_get(m_work, r - 1, c)):
                if inject_seq is not None:
                    _set(m_work, r - 1, c, inject_seq[num])
                num += 1
                _set(n, r - 1, c, num)
                seq.append(_get(m_work, r - 1, c))
            r -= 1

        elif d == 'down' and VERT_MOVE.match(_get(m_work, r + 1, c)):
            if NT_MATCH.match(_get(m_work, r + 1, c)):
                if inject_seq is not None:
                    _set(m_work, r + 1, c, inject_seq[num])
                num += 1
                _set(n, r + 1, c, num)
                seq.append(_get(m_work, r + 1, c))
            r += 1

        # ── Corner crossovers ─────────────────────────────────────────
        d, r, c = _apply_crossover(m_work, d, r, c)

        # ── Check for 3′ end ──────────────────────────────────────────
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nb = _get(m_work, r + dr, c + dc)
            if DIGIT.match(nb) and int(nb) == 3:
                success = True
                break
        if success:
            break

    if not success:
        print(f"The trace through the structure failed (3′ end not found). "
              f"Ended at row {r}, column {c}.")

    return seq, n, r, c, d, success


# ── Base-pair finder ──────────────────────────────────────────────────────────

def find_base_pairs(
    m: Grid, n: Grid, r0: int, c0: int, d0: str, nt: int,
    kl_pat: List[str],
    strand_dir: Optional[Grid] = None,
) -> Tuple[List[int], List[int], List[str]]:
    """
    Second pass through the strand: record pairing relationships.

    p-values encode the pair type:
        'p' or 'b' → standard Watson-Crick (→ '(' / ')')
        '!'        → forced/extra pair (→ '{' / '}')
        '*'        → kissing-loop pair  (→ '[' / ']')
        KL char    → from @-pattern     (printed as-is)
        '-' or 'i' → unpaired           (→ '.')

    Returns parallel arrays (a, b, p) where a[i] and b[i] are 1-based nt
    positions, and p[i] is the pair type.  b[i]==0 means unpaired.
    """
    a_arr: List[int] = []
    b_arr: List[int] = []
    p_arr: List[str] = []

    r, c = r0, c0
    d    = d0
    current_kl = 0

    for _ in range(nt + 10000):
        # ── Horizontal strand (right / left) ──────────────────────────
        if d == 'right' and HORIZ_MOVE.match(_get(m, r, c + 1)):
            if NT_MATCH.match(_get(m, r, c + 1)):
                _record_pair_horiz(m, n, a_arr, b_arr, p_arr, kl_pat,
                                   r, c + 1, d, 1, strand_dir, current_kl)
                if _get(m, r, c + 1) == 'X' and not _get(m, r, c + 2) == 'X':
                    current_kl += 1
            c += 1

        elif d == 'left' and HORIZ_MOVE.match(_get(m, r, c - 1)):
            if NT_MATCH.match(_get(m, r, c - 1)):
                _record_pair_horiz(m, n, a_arr, b_arr, p_arr, kl_pat,
                                   r, c - 1, d, -1, strand_dir, current_kl)
                if _get(m, r, c - 1) == 'X' and not _get(m, r, c - 2) == 'X':
                    current_kl += 1
            c -= 1

        # ── Vertical strand (down / up) ───────────────────────────────
        elif d == 'down' and VERT_MOVE.match(_get(m, r + 1, c)):
            if NT_MATCH.match(_get(m, r + 1, c)):
                _record_pair_vert(m, n, a_arr, b_arr, p_arr,
                                  r + 1, c, d, strand_dir)
            r += 1

        elif d == 'up' and VERT_MOVE.match(_get(m, r - 1, c)):
            if NT_MATCH.match(_get(m, r - 1, c)):
                _record_pair_vert(m, n, a_arr, b_arr, p_arr,
                                  r - 1, c, d, strand_dir)
            r -= 1

        # ── Blank horizontal / vertical connectors ────────────────────
        if d == 'right' and HORIZ_STEP.match(_get(m, r, c + 1)):
            c += 1
        elif d == 'left' and HORIZ_STEP.match(_get(m, r, c - 1)):
            c -= 1
        elif d == 'down' and VERT_STEP.match(_get(m, r + 1, c)):
            r += 1
        elif d == 'up' and VERT_STEP.match(_get(m, r - 1, c)):
            r -= 1

        # ── Corner crossovers ─────────────────────────────────────────
        d, r, c = _apply_crossover(m, d, r, c)

        # ── Check for 3′ end ──────────────────────────────────────────
        done = False
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nb = _get(m, r + dr, c + dc)
            if DIGIT.match(nb) and int(nb) == 3:
                done = True
                break
        if done:
            break

    return a_arr, b_arr, p_arr


def _record_pair_horiz(m, n, a_arr, b_arr, p_arr, kl_pat,
                        r, c, d, step, strand_dir, current_kl):
    """Record a pairing entry for a horizontal strand cell."""
    below = _get(m, r + 1, c)
    above = _get(m, r - 1, c)
    if PAIR_MARKER.match(below):
        a_arr.append(_get(n, r, c) or 0)
        b_arr.append(_get(n, r + 2, c) or 0)
        p_arr.append(below)
        if strand_dir is not None:
            _set(strand_dir, r + 1, c, d)
    elif PAIR_MARKER.match(above):
        a_arr.append(_get(n, r, c) or 0)
        b_arr.append(_get(n, r - 2, c) or 0)
        p_arr.append(above)
        if strand_dir is not None:
            _set(strand_dir, r - 1, c, d)
    else:
        a_arr.append(_get(n, r, c) or 0)
        if _get(m, r, c) == 'X' and kl_pat and current_kl < len(kl_pat):
            b_arr.append(0)
            p_arr.append(kl_pat[current_kl])
        else:
            b_arr.append(0)
            p_arr.append('-')


def _record_pair_vert(m, n, a_arr, b_arr, p_arr, r, c, d, strand_dir):
    """Record a pairing entry for a vertical strand cell."""
    right = _get(m, r, c + 1)
    left  = _get(m, r, c - 1)
    if BPAIR_MARKER.match(right):
        a_arr.append(_get(n, r, c) or 0)
        b_arr.append(_get(n, r, c + 2) or 0)
        p_arr.append(right)
        if strand_dir is not None:
            _set(strand_dir, r, c + 1, d)
    elif BPAIR_MARKER.match(left):
        a_arr.append(_get(n, r, c) or 0)
        b_arr.append(_get(n, r, c - 2) or 0)
        p_arr.append(left)
        if strand_dir is not None:
            _set(strand_dir, r, c - 1, d)
    else:
        a_arr.append(_get(n, r, c) or 0)
        b_arr.append(0)
        p_arr.append('i')


# ── Crossover helper ──────────────────────────────────────────────────────────

def _apply_crossover(m: Grid, d: str, r: int, c: int) -> Tuple[str, int, int]:
    """Apply a corner turn if the next cell is a crossover symbol."""
    if d == 'right':
        nxt = _get(m, r, c + 1)
        if nxt == '7': return 'down', r, c + 1
        if nxt == 'J': return 'up',   r, c + 1
    elif d == 'left':
        nxt = _get(m, r, c - 1)
        if nxt == 'L': return 'up',    r, c - 1
        if nxt == 'r': return 'down',  r, c - 1
    elif d == 'down':
        nxt = _get(m, r + 1, c)
        if nxt == 'L': return 'right', r + 1, c
        if nxt == 'J': return 'left',  r + 1, c
    elif d == 'up':
        nxt = _get(m, r - 1, c)
        if nxt == '7': return 'left',  r - 1, c
        if nxt == 'r': return 'right', r - 1, c
    return d, r, c


# ── Structure string builder ──────────────────────────────────────────────────

def build_structure_string(a_arr: List[int], b_arr: List[int],
                            p_arr: List[str]) -> Tuple[str, str]:
    """
    Convert (a, b, p) arrays into a dot-bracket string and a simplified
    structure_map string (PKs treated as plain pairs/dots).

    Returns (dot_bracket, structure_map).
    """
    dot_bracket  = []
    structure_map = []

    _OPEN  = {'p': '(', 'b': '(', '!': '{', '*': '['}
    _CLOSE = {'p': ')', 'b': ')', '!': '}', '*': ']'}
    _OPEN_MAP  = {'p': '(', 'b': '(', '!': '(', '*': '['}
    _CLOSE_MAP = {'p': ')', 'b': ')', '!': ')', '*': ']'}

    for idx, ai in enumerate(a_arr):
        bi = b_arr[idx]
        pi = p_arr[idx]

        if bi != 0:
            if ai > bi:          # closing bracket
                ch    = _CLOSE.get(pi, '.')
                ch_sm = _CLOSE_MAP.get(pi, ')')
            else:                # opening bracket
                ch    = _OPEN.get(pi, '.')
                ch_sm = _OPEN_MAP.get(pi, '(')
        else:
            ch    = pi if KL_CHARS.match(pi) else '.'
            ch_sm = '.'

        dot_bracket.append(ch)
        structure_map.append(ch_sm)

    return ''.join(dot_bracket), ''.join(structure_map)


# ── Grid renderer ─────────────────────────────────────────────────────────────

def render_grid_base(m: Grid, strand_dir: Optional[Grid], cross: Grid,
                     show_nts: bool = True) -> str:
    """
    Render the grid to a string with Unicode box-drawing characters.
    Nucleotide cells are shown as their symbol if show_nts=True, else as '─' or '│'.
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

            if ch == '\n':
                out.append('\n')
                continue
            if ch == ' ':
                out.append(' ')
                continue
            if ch in ('3', '5'):
                out.append(ch)
                continue
            if ch == 'T':
                out.append('─')
                continue

            # Box-drawing symbols
            sym = SYMBOL_TO_UNICODE.get(ch)
            if sym:
                out.append(sym)
                continue

            # Pair connectors — direction-dependent
            if ch in ('p', '!', '*'):
                sd = _get(strand_dir, i, j) if strand_dir else None
                if ch == 'p':
                    out.append('┊')
                elif ch == '!':
                    out.append('!')
                elif ch == '*':
                    out.append('*')
                continue

            if ch == 'b':
                out.append('─')
                continue

            # Cross
            if cross and _get(cross, i, j) == '^':
                out.append('^')
                continue
            if ch == 'i':
                out.append('│')
                continue
            if ch == 'x':
                out.append('┼')
                continue

            # Nucleotide cells
            if NT_MATCH.match(ch):
                if show_nts:
                    out.append(ch)
                else:
                    sd = _get(strand_dir, i, j) if strand_dir else None
                    if sd in ('down', 'up'):
                        out.append('│')
                    else:
                        out.append('─')
                continue

            out.append(ch)

    return ''.join(out)


# ── Pair-map builder (for sequence scrubbing in trace_analysis) ────────────────

def build_pair_map_from_structure(structure: str) -> List[int]:
    """Build an array where pair_map[i] = j means i pairs with j (0-indexed)."""
    n = len(structure)
    pair_map = list(range(n))
    bracket_stack: Dict[str, List[int]] = {'(': [], '[': []}

    for idx, ch in enumerate(structure):
        if ch == '(':
            bracket_stack['('].append(idx)
        elif ch == ')':
            if bracket_stack['(']:
                partner = bracket_stack['('].pop()
                pair_map[idx] = partner
                pair_map[partner] = idx
        elif ch == '[':
            bracket_stack['['].append(idx)
        elif ch == ']':
            if bracket_stack['[']:
                partner = bracket_stack['['].pop()
                pair_map[idx] = partner
                pair_map[partner] = idx

    return pair_map
