#!/usr/bin/env python3
"""
Dragon v2.03 - RNA Sequence Designer (Python Port)
Original Perl: Copyright 2020 Cody Geary and Ebbe S. Andersen
Python port: translated from dragon.pl

Port of dragon.pl to Python.

Behavior note:
  This file now preserves dragon.pl control flow, including the original
  drift-counter semantics, so the mutation cycle matches the Perl program
  as closely as possible while still using ViennaRNA Python bindings.

Usage:
    python dragon.py [output_directory]

Reads target.txt from the current directory.
Requires ViennaRNA Python bindings: pip install ViennaRNA
"""

import sys
import os
import random
import time
import subprocess
from typing import List

try:
    import RNA
except ImportError:
    sys.exit(
        "ViennaRNA Python bindings not found.\n"
        "Install with:  pip install ViennaRNA\n"
        "Or see https://www.tbi.univie.ac.at/RNA/"
    )

# ============================================================
# User-Configurable Parameters
# ============================================================

INIT_SETUP          = 5      # Cycles of the initial setup phase
SITE_MUTATION_RATE  = 45     # Rate of mutation in mispaired regions
MUTATE_LOOP         = 2      # 1-100, rate at which mispaired loop nts are mutated
DRIFT_PERIOD        = 3      # Fire a drift event every Nth round
MY_DRIFT_RATE       = 2      # Rate for random point-mutations to the entire sequence
GC_MUT_RATIO        = 15     # How much GCs are preferred over AUs
DOUBLE_MUTANT       = 99     # Percent of mutants targeting base-pairs vs single positions
MISFOLD_TOLERANCE   = 1.06   # Factor of allowed distance increase during mutation
OUTPUT_THRESHOLD    = 0      # 1=verbose, 0=hide output for neutral mutations
MAX_GC              = 60     # Desired GC content of initial design
MAX_GC_OPT          = 58     # Final GC content target after optimization
KL_GC               = 30     # GC% for Kissing Loops
TARGET_GU           = 1      # Minimum GU wobble pairs
LONGRANGE           = 250    # Min distance for pairs requiring AUC-limited alphabet

DUPLICATE_WINDOW    = 9      # Length of duplicate regions Pattern Search looks for
COMPLEMENT_WINDOW   = 10     # Window size to search for complement sequences
PALINDROME_WINDOW   = 8      # Min palindrome size to search and mark

MAX_FAILED_TRIES    = 100    # Times to try mutating surrounding nts when locked by constraints

GENERATE_PREVIEWS   = False  # Run trace_analysis preview every new round

KL_MIN_DELTA_G      = -6.0   # Minimum threshold for non-cognate KL interactions
MIN_KL              = -9.2   # KL scoring: below this energy is too weak
MAX_KL              = -11.0  # KL scoring: above this energy is too strong

FAV_RAD_LEVEL       = 15     # ~Number of mutants per round (adaptive)
TARGET_ED           = 1      # Target ensemble diversity to stop optimization

# ============================================================
# ViennaRNA API Wrappers
# ============================================================

def rna_fold(sequence: str) -> str:
    """Return MFE dot-bracket structure for a sequence."""
    structure, _mfe = RNA.fold(sequence)
    return structure


def rna_distance(fold1: str, fold2: str) -> int:
    """Return base-pair distance between two dot-bracket structures."""
    return RNA.bp_distance(fold1.strip(), fold2.strip())


def rna_duplex(seq1: str, seq2: str) -> tuple:
    """
    Return (dot-bracket structure, free energy) for duplex folding of two sequences.
    Replaces: qx(rnaduplex < KLtest.txt) and the broken qr// regex parse.
    """
    duplex = RNA.duplexfold(seq1, seq2)
    return duplex.structure, float(duplex.energy)


def rna_ensemble_diversity(sequence: str) -> float:
    """
    Compute ensemble diversity via partition function.
    Replaces: qx(rnafold --noPS -p < seq.txt) and the broken qr// parse of its output.
    """
    fc = RNA.fold_compound(sequence)
    _structure, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()
    return fc.mean_bp_distance()

# ============================================================
# Sequence Complement Utilities
# ============================================================

_COMPLEMENT = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}


def reverse_complement(seq: str) -> str:
    return ''.join(_COMPLEMENT.get(b, 'N') for b in reversed(seq))


def rand_base() -> str:
    return random.choice(['A', 'U', 'C', 'G'])


def rand_purine() -> str:
    return random.choice(['A', 'G'])


def rand_pyrimidine() -> str:
    return random.choice(['C', 'U'])


def rand_strong() -> str:
    return random.choice(['C', 'G'])


def rand_weak() -> str:
    return random.choice(['A', 'U'])

# ============================================================
# Pair Map Builder
# (renamed from 'map' in Perl, which shadowed the built-in)
# ============================================================

def build_pair_map(secondary: str) -> List[int]:
    """
    Build a pairing map from a dot-bracket structure, including pseudoknots.
    Returns a list where pair_map[i] = j means i pairs with j.
    Unpaired positions map to themselves: pair_map[i] = i.
    Handles: ( ) [ ] { } A-J (open) / 1-0 (close)
    """
    n = len(secondary)
    pair_map = list(range(n))

    open_to_close = {
        '(': ')', '[': ']', '{': '}',
        'A': '1', 'B': '2', 'C': '3', 'D': '4', 'E': '5',
        'F': '6', 'G': '7', 'H': '8', 'I': '9', 'J': '0',
    }
    close_to_open = {v: k for k, v in open_to_close.items()}
    stacks = {op: [] for op in open_to_close}

    for idx, ch in enumerate(secondary):
        if ch == '.':
            pair_map[idx] = idx
        elif ch in open_to_close:
            stacks[ch].append(idx)
        elif ch in close_to_open:
            opener = close_to_open[ch]
            if stacks[opener]:
                partner = stacks[opener].pop()
                pair_map[idx] = partner
                pair_map[partner] = idx

    return pair_map

# ============================================================
# Input Parsing
# ============================================================

def read_target(filename: str = 'target.txt'):
    """
    Read target.txt. Returns (name, target_pk_str, constraint_str, start_seq_str).
    Line 0: design name
    Line 1: target structure with pseudoknot notation
    Line 2: constraint string
    Line 3: (optional) starting sequence
    """
    with open(filename, 'r') as f:
        content = f.read()
    content = content.replace('\r\n', '\n').replace('\r', '\n')
    lines = content.split('\n')

    name          = lines[0].strip() if len(lines) > 0 else 'design'
    target_pk_str = lines[1].strip() if len(lines) > 1 else ''
    constraint_str= lines[2].strip() if len(lines) > 2 else ''
    start_seq_str = lines[3].strip() if len(lines) > 3 else ''
    return name, target_pk_str, constraint_str, start_seq_str


# Characters that denote pseudoknot structure (not standard Watson-Crick)
_PK_CHARS      = set('[]{}ABCDEFGHIJ1234567890')
_OPEN_PK       = set('{ABCDEFGHIJ')
_CLOSE_PK      = set('}1234567890')


def parse_target(target_pk_str: str):
    """
    Derive three representations of the target structure:
      init_target    : original character list (preserved for '{' dovetail detection)
      target_nopk    : all PK chars replaced with '.' (used for distance calculation)
      target_scrubbed: exotic PK brackets unified to [ ] (used for pair_map and main logic)
    """
    init_target = list(target_pk_str)
    target_nopk_chars = []
    target_scrubbed_chars = []

    for ch in target_pk_str:
        target_nopk_chars.append('.' if ch in _PK_CHARS else ch)
        if ch in _OPEN_PK:
            target_scrubbed_chars.append('[')
        elif ch in _CLOSE_PK:
            target_scrubbed_chars.append(']')
        else:
            target_scrubbed_chars.append(ch)

    return init_target, list(target_nopk_chars), list(target_scrubbed_chars)

# ============================================================
# GC Content Counter
# ============================================================

def count_gc(sol: List[str], target: List[str], init_target: List[str],
             pair_map: List[int], strand_length: int) -> dict:
    """Count pair types and compute GC%, KL GC%. Returns a stats dict."""
    pairs = gc_pairs = au_pairs = gu_pairs = 0
    kl_pairs = kl_gc_pairs = 0

    for i in range(strand_length):
        if target[i] == '(':
            pairs += 1
            a, b = sol[i], sol[pair_map[i]]
            if {a, b} == {'C', 'G'}:
                gc_pairs += 1
            elif {a, b} == {'A', 'U'}:
                au_pairs += 1
            elif {a, b} == {'G', 'U'}:
                gu_pairs += 1
        elif init_target[i] == '[':
            kl_pairs += 1
            a, b = sol[i], sol[pair_map[i]]
            if {a, b} == {'C', 'G'}:
                kl_gc_pairs += 1

    gc_cont    = round(1000 * gc_pairs / pairs) / 10 if pairs > 0 else -1
    kl_gc_cont = round(1000 * kl_gc_pairs / kl_pairs) / 10 if kl_pairs > 0 else 0

    return dict(pairs=pairs, gc_pairs=gc_pairs, au_pairs=au_pairs,
                gu_pairs=gu_pairs, kl_pairs=kl_pairs, kl_gc_pairs=kl_gc_pairs,
                gc_cont=gc_cont, kl_gc_cont=kl_gc_cont)

# ============================================================
# Pattern / Repeat Detector
# ============================================================

# Restriction sites to avoid (RNA notation, U not T)
_RESTRICTION_SITES = [
    'GGUCUC',   # BsaI
    'GAGACC',   # BsaI
    'GAAGAC',   # BbsI
    'GUCUUC',   # BbsI
    'CGUCUC',   # BsMBI
    'GAGACG',   # BsMBI
    'UUUU',     # Pol III terminator
    'GCUCUUC',  # SapI
    'GAAGAGC',  # SapI
    'AUCUGUU',  # PTH
]


def count_repeats(trial_sol: List[str], target: List[str], constraint: List[str],
                  strand_length: int,
                  complement_window: int, duplicate_window: int,
                  palindrome_window: int) -> dict:
    """
    Scan for sequence patterns that should be avoided.
    repeat_map key:
      '-'  clean
      'P'  complementary region (potential off-target duplex)
      'L'  palindrome
      'D'  duplicated sequence
      'W'  AU-rich run of 6+
      'S'  GC-rich run of 6+
      'G'/'C'/'A'/'U'  poly-homopolymer run of 5+
      'X'  restriction enzyme site
    """
    seq = ''.join(trial_sol)
    repeat_map = ['-'] * strand_length

    # --- Complementary regions ---
    for i in range(strand_length - complement_window):
        window = seq[i:i + complement_window]
        rc = reverse_complement(window)
        start = 0
        while True:
            pos = seq.find(rc, start)
            if pos == -1:
                break
            for k in range(complement_window):
                if i + k < strand_length:
                    repeat_map[i + k] = 'P'
                if pos + k < strand_length:
                    repeat_map[pos + k] = 'P'
            start = pos + 1

    complement_zones = sum(1 for ch in repeat_map if ch == 'P')

    # --- Palindromes ---
    for i in range(strand_length - palindrome_window):
        window = seq[i:i + palindrome_window]
        rc = reverse_complement(window)
        start = i
        while start < strand_length:
            pos_rc = seq.find(rc, start)
            pos_s  = seq.find(window, start)
            if pos_rc == -1 or pos_s == -1:
                break
            if pos_rc == pos_s:  # sequence is its own reverse complement: palindrome
                for k in range(palindrome_window):
                    if i + k < strand_length:
                        repeat_map[i + k] = 'L'
                    if pos_rc + k < strand_length:
                        repeat_map[pos_rc + k] = 'L'
            start = max(pos_rc, pos_s) + 1

    palindromes = sum(1 for ch in repeat_map if ch == 'L')

    # --- Duplicated sequences ---
    for i in range(strand_length - duplicate_window):
        window = seq[i:i + duplicate_window]
        start = i
        while True:
            pos = seq.find(window, start)
            if pos == -1:
                break
            if pos > i:
                for k in range(duplicate_window):
                    if i + k < strand_length:
                        repeat_map[i + k] = 'D'
                    if pos + k < strand_length:
                        repeat_map[pos + k] = 'D'
            start = pos + 1

    duplication_zones = sum(1 for ch in repeat_map if ch == 'D')

    # --- AU-rich and GC-rich runs of 6+ ---
    pattern_repeats = 0
    for i in range(strand_length - 5):
        window6 = trial_sol[i:i + 6]
        if all(b in ('A', 'U') for b in window6):
            pattern_repeats += 1
            for k in range(6):
                repeat_map[i + k] = 'W'
        if all(b in ('C', 'G') for b in window6):
            if all(constraint[i + k] == 'S' for k in range(6)):
                print(f"\nLocked constraint of 6S found at position {i}")
            else:
                pattern_repeats += 1
            for k in range(6):
                repeat_map[i + k] = 'S'

    # --- Poly-homopolymer runs of 5+ ---
    poly_repeats = 0
    for i in range(strand_length - 6):
        if all(trial_sol[i + k] == trial_sol[i] for k in range(1, 6)):
            poly_repeats += 1
            base = trial_sol[i]
            for k in range(6):
                if i + k < strand_length:
                    repeat_map[i + k] = base

    # --- Restriction sites ---
    restriction_sites = 0
    for i in range(strand_length - 3):
        for site in _RESTRICTION_SITES:
            end = i + len(site)
            if end <= strand_length and seq[i:end] == site:
                restriction_sites += 1
                for k in range(len(site)):
                    repeat_map[i + k] = 'X'
                break

    return dict(
        repeat_map=repeat_map,
        complement_zones=complement_zones,
        palindromes=palindromes,
        duplication_zones=duplication_zones,
        pattern_repeats=pattern_repeats,
        poly_repeats=poly_repeats,
        restriction_sites=restriction_sites,
    )

# ============================================================
# KL Repeat Detector
# ============================================================

def pattern_prevent(trial_sol: List[str], target_nopk: List[str],
                    constraint: List[str], strand_length: int) -> tuple:
    """
    Search for repeated kissing-loop sequences (180KL and 120KL patterns).
    Returns (pattern_zones, kl_repeats).
    """
    pattern_zones = ['-'] * strand_length
    kl_repeats = 0

    for i in range(strand_length - 10):
        for j in range(i + 1, strand_length - 9):

            # 180KL: 8-nt match in two unpaired regions
            if (all(trial_sol[i + k] == trial_sol[j + k] for k in range(8)) and
                    all(target_nopk[i + k] == '.' for k in range(8)) and
                    all(target_nopk[j + k] == '.' for k in range(8)) and
                    constraint[j + 4] == 'N'):
                for k in range(8):
                    pattern_zones[i + k] = trial_sol[i + k]
                    pattern_zones[j + k] = trial_sol[i + k]
                kl_repeats += 1

            # 120KL: (.......)  loop, 7-nt internal match
            elif (j + 8 < strand_length and
                  all(trial_sol[i + k + 1] == trial_sol[j + k + 1] for k in range(7)) and
                  target_nopk[i] == '(' and
                  all(target_nopk[i + k + 1] == '.' for k in range(7)) and
                  target_nopk[i + 8] == ')' and
                  target_nopk[j] == '(' and
                  all(target_nopk[j + k + 1] == '.' for k in range(7)) and
                  target_nopk[j + 8] == ')' and
                  constraint[j + 4] == 'N'):
                for k in range(1, 8):
                    pattern_zones[i + k] = trial_sol[i + k]
                    pattern_zones[j + k] = trial_sol[i + k]
                kl_repeats += 1

    return pattern_zones, kl_repeats

# ============================================================
# KL Thermodynamic Analyser
# ============================================================

def analyze_kl(trial_sol: List[str], target: List[str], target_nopk: List[str],
               pair_map: List[int], constraint: List[str], strand_length: int,
               duplex_cache: dict, report_kl: bool, full_kl_output: bool,
               kl_min_delta_g: float, min_kl: float, max_kl: float,
               log) -> dict:
    """
    Find all kissing loops, score their thermodynamics, check non-cognate pairings.
    Returns a results dict with kl_score, lists of sequences/energies, opt targets, etc.
    """
    num_kl = 0
    kl_score = 0
    kl_notcounted = 0
    # These are 1-indexed arrays (slot 0 unused) matching the Perl convention
    a_list       = [None]
    b_list       = [None]
    e_list       = [None]
    nt_list      = [None]
    kl_type_list = [None]
    kl_specified = [None]
    kl_opt_target = [0] * strand_length
    num_kl_sites  = 0
    non_cognates  = ''

    def do_duplex(s1: str, s2: str) -> tuple:
        """Memoised rna_duplex call."""
        key = f"{s1}.{s2}"
        key_r = f"{s2}.{s1}"
        if key in duplex_cache:
            return duplex_cache[key]
        if key_r in duplex_cache:
            return duplex_cache[key_r]
        result = rna_duplex(s1, s2)
        duplex_cache[key] = result
        return result

    if report_kl:
        log("Mapping Kissing Loops:\n")

    # --- Identify KL sites ---
    for i in range(strand_length - 8):

        # 180KL: .]]]]]]  or  )]]]]]]  followed by .
        if ((target[i] in ('.', ')')) and
                target[i+1:i+7] == [']']*6 and
                target[i+7] == '.'):

            bkl_flag = 1 if target[i] == ')' else 0
            kl_b = ''.join(trial_sol[i + (1 if bkl_flag else 0): i + 8])

            m1 = pair_map[i + 1]
            if m1 > 0 and target[m1 - 1] == '(':
                bkl_flag = 2
                m6 = pair_map[i + 6]
                kl_a = ''.join([
                    trial_sol[m6 - 1], trial_sol[m6],
                    trial_sol[pair_map[i+5]], trial_sol[pair_map[i+4]],
                    trial_sol[pair_map[i+3]], trial_sol[pair_map[i+2]],
                    trial_sol[pair_map[i+1]],
                ])
            else:
                m6 = pair_map[i + 6]
                kl_a = ''.join([
                    trial_sol[m6 - 1], trial_sol[m6],
                    trial_sol[pair_map[i+5]], trial_sol[pair_map[i+4]],
                    trial_sol[pair_map[i+3]], trial_sol[pair_map[i+2]],
                    trial_sol[pair_map[i+1]], trial_sol[pair_map[i+1] + 1],
                ])

            num_kl += 1
            struct, energy = do_duplex(kl_a, kl_b)

            if report_kl:
                label = {0: '180 KL', 1: 'KL-BKL', 2: 'BKL-KL'}.get(bkl_flag, '180 KL')
                log(f"{'  ' if num_kl < 10 else ''}{num_kl},{label},"
                    f"{kl_a},{kl_b},{struct},{energy}")

            a_list.append(kl_a)
            b_list.append(kl_b)
            e_list.append(energy)
            nt_list.append(i)
            kl_type_list.append('A')

            all_n = all(constraint[i + k] == 'N' for k in range(1, 7))
            if all_n:
                kl_specified.append('N')
            else:
                kl_specified.append('Y')
                kl_notcounted += 1
                if report_kl:
                    log(" (ignored)")

            if report_kl:
                log("\n")

            # Score: outside the desired energy window → penalty
            if (energy > min_kl or energy < max_kl) and kl_specified[-1] == 'N':
                kl_score += 2
                num_kl_sites += 6
                for k in range(1, 7):
                    if i + k < strand_length:
                        kl_opt_target[i + k] += 1

        # 120KL: [[[[[[[
        elif (i + 6 < strand_length and
              all(target[i + k] == '[' for k in range(7))):

            kl_a = ''.join(trial_sol[i:i + 7])
            kl_b = ''.join([
                trial_sol[pair_map[i+6]], trial_sol[pair_map[i+5]],
                trial_sol[pair_map[i+4]], trial_sol[pair_map[i+3]],
                trial_sol[pair_map[i+2]], trial_sol[pair_map[i+1]],
                trial_sol[pair_map[i]],
            ])

            num_kl += 1
            struct, energy = do_duplex(kl_a, kl_b)

            if report_kl:
                log(f"{'  ' if num_kl < 10 else ''}{num_kl},120 KL,"
                    f"{kl_a},{kl_b},{struct},{energy}")

            a_list.append(kl_a)
            b_list.append(kl_b)
            e_list.append(energy)
            nt_list.append(i)
            kl_type_list.append('B')

            all_n = all(constraint[i + k] == 'N' for k in range(7))
            if all_n:
                kl_specified.append('N')
            else:
                kl_specified.append('Y')
                kl_notcounted += 1
                if report_kl:
                    log(" (ignored)")

            if report_kl:
                log("\n")

            if (energy > min_kl or energy < max_kl) and kl_specified[-1] == 'N':
                kl_score += 2
                num_kl_sites += 7
                for k in range(7):
                    if i + k < strand_length:
                        kl_opt_target[i + k] += 1

    calculation = kl_score / 2
    log(f"\n\n{calculation} KLs with energy lower than {min_kl} kcal "
        f"or greater than {max_kl} kcal were counted.\n\n")
    if kl_notcounted > 0:
        log(f"{kl_notcounted} KLs not counted (user-specified).\n\n")

    # --- Non-cognate pairings ---
    header = (f"\n\nMapping Non-Cognate Pairings (non-cognate at least as strong as "
              f"on-target, or stronger than {kl_min_delta_g} kcal):\n")
    if report_kl:
        log(header)
    if full_kl_output:
        log(f"\n\nMapping Non-Cognate Pairings (all non-cognate at least as strong "
            f"as {kl_min_delta_g} kcal):\n")
    non_cognates = header

    def _add_targets(idx: int):
        """Add positions to the optimization target list for KL idx."""
        nonlocal num_kl_sites
        ktype = kl_type_list[idx]
        nt = nt_list[idx]
        if ktype == 'A':
            num_kl_sites += 6
            for k in range(1, 7):
                if nt + k < strand_length:
                    kl_opt_target[nt + k] += 1
        elif ktype == 'B':
            num_kl_sites += 7
            for k in range(7):
                if nt + k < strand_length:
                    kl_opt_target[nt + k] += 1

    for i in range(1, num_kl):
        for j in range(i + 1, num_kl + 1):
            for s1, s2, label in [
                (a_list[i], a_list[j], f"{i}(A) + {j}(A)"),
                (a_list[i], b_list[j], f"{i}(A) + {j}(B)"),
                (b_list[i], a_list[j], f"{i}(B) + {j}(A)"),
                (b_list[i], b_list[j], f"{i}(B) + {j}(B)"),
            ]:
                struct, energy = do_duplex(s1, s2)
                if energy <= e_list[i] or energy <= e_list[j] or energy <= kl_min_delta_g:
                    if kl_specified[i] == 'N' or kl_specified[j] == 'N':
                        kl_score += 1
                    _add_targets(i)
                    _add_targets(j)
                should_print = (
                    (full_kl_output and energy < kl_min_delta_g) or
                    (report_kl and (
                        energy <= e_list[i] or
                        energy <= e_list[j] or
                        energy <= kl_min_delta_g
                    ))
                )
                if should_print:
                    line = f"{label},{s1},{s2},{struct},{energy}\n"
                    log(line)
                    non_cognates += line

    log(f"\nFinal KL Score = {kl_score}\n")

    return dict(
        num_kl=num_kl, kl_score=kl_score, kl_notcounted=kl_notcounted,
        a_list=a_list, b_list=b_list, e_list=e_list,
        nt_list=nt_list, kl_type=kl_type_list, kl_specified=kl_specified,
        kl_opt_target=kl_opt_target, num_kl_sites=num_kl_sites,
        non_cognates=non_cognates,
    )

# ============================================================
# Sequence Initialiser
# ============================================================

def initialize_sequence(target: List[str], init_target: List[str],
                        constraint: List[str], pair_map: List[int],
                        strand_length: int) -> List[str]:
    """Build an initial sequence satisfying constraints and target structure."""
    sol = ['X'] * strand_length

    # --- Unpaired (loop) positions ---
    for i in range(strand_length):
        if target[i] == '.':
            c = constraint[i]
            if c in ('A', 'U', 'C', 'G'):
                sol[i] = c
            elif c == 'R':
                sol[i] = rand_purine()
            elif c == 'N':
                sol[i] = rand_base()
            elif c == 'Y':
                sol[i] = rand_pyrimidine()
            elif c == 'S':
                sol[i] = rand_strong()
            elif c == 'W':
                sol[i] = rand_weak()
            else:
                sol[i] = 'A'

    # --- Stem positions ---
    for i in range(strand_length):
        if target[i] not in ('(', '['):
            continue
        ci  = constraint[i]
        mi  = pair_map[i]
        cm  = constraint[mi]

        if ci == 'A':
            sol[i] = 'A'; sol[mi] = 'U'
        elif ci == 'C':
            sol[i] = 'C'; sol[mi] = 'G'
        elif ci == 'U':
            sol[i] = 'U'
            sol[mi] = 'A' if cm in ('A', 'N') else 'G' if cm in ('G', 'K') else sol[mi]
        elif ci == 'G':
            sol[i] = 'G'
            sol[mi] = 'C' if cm in ('C', 'N', 'Y') else 'U' if cm in ('U', 'K') else sol[mi]
        elif ci == 'K':
            if cm in ('G', 'R', 'N'):
                sol[mi] = 'U'
            elif cm == 'U':
                sol[mi] = 'G'
            if cm == 'K':
                bp_dir = random.randint(0, 1)
                if (mi - i) > LONGRANGE:
                    bp_dir = 1
                if bp_dir == 0 or sol[i] == 'G' or sol[mi] == 'U':
                    sol[i] = 'G'; sol[mi] = 'U'
                else:
                    sol[i] = 'U'; sol[mi] = 'G'
        elif ci == 'S':
            if cm == 'G':
                sol[mi] = 'C'
            elif cm == 'C':
                sol[mi] = 'G'
            if cm == 'S':
                bp_dir = random.randint(0, 1)
                if (mi - i) > LONGRANGE:
                    bp_dir = 1
                if bp_dir == 0 or sol[i] == 'G' or sol[mi] == 'C':
                    sol[i] = 'G'; sol[mi] = 'C'
                else:
                    sol[i] = 'C'; sol[mi] = 'G'
        elif ci == 'R':
            if cm == 'U':
                sol[mi] = 'A'
            elif cm == 'C':
                sol[mi] = 'G'
            if cm in ('Y', 'N'):
                bp_dir = random.randint(0, 1)
                if bp_dir == 0 or sol[i] == 'G' or sol[mi] == 'C':
                    sol[i] = 'G'; sol[mi] = 'C'
                else:
                    sol[i] = 'A'; sol[mi] = 'U'
        elif ci == 'Y':
            if cm == 'G':
                sol[mi] = 'C'
            elif cm == 'A':
                sol[mi] = 'U'
            if cm in ('R', 'N'):
                bp_dir = random.randint(0, 1)
                if bp_dir == 0 or sol[i] == 'U' or sol[mi] == 'A':
                    sol[i] = 'C'; sol[mi] = 'G'
                else:
                    sol[i] = 'U'; sol[mi] = 'A'

        # N:N standard stems
        if target[i] == '(' and ci == 'N':
            if cm == 'A':
                sol[i] = 'U'; sol[mi] = 'A'
            elif cm == 'U':
                sol[i] = 'A'; sol[mi] = 'U'
            elif cm == 'C':
                sol[i] = 'G'; sol[mi] = 'C'
            elif cm == 'G':
                sol[i] = 'C'; sol[mi] = 'G'
            elif cm == 'N' and sol[i] == 'X' and sol[mi] == 'X':
                bp_dir  = random.randint(0, 1)
                bp_type = random.uniform(0, 100)
                if bp_type < (100 - MAX_GC - TARGET_GU / 2):
                    sol[i], sol[mi] = ('A', 'U') if bp_dir == 0 else ('U', 'A')
                else:
                    sol[i], sol[mi] = ('G', 'C') if bp_dir == 0 else ('C', 'G')

        # N:N KL stems
        if target[i] == '[' and ci == 'N':
            if cm == 'A':
                sol[i] = 'U'; sol[mi] = 'A'
            elif cm == 'U':
                sol[i] = 'A'; sol[mi] = 'U'
            elif cm == 'C':
                sol[i] = 'G'; sol[mi] = 'C'
            elif cm == 'G':
                sol[i] = 'C'; sol[mi] = 'G'
            elif cm == 'N' and sol[i] == 'X' and sol[mi] == 'X':
                bp_dir  = random.randint(0, 1)
                bp_type = random.uniform(0, 100)
                if bp_type < KL_GC:
                    sol[i], sol[mi] = ('G', 'C') if bp_dir == 0 else ('C', 'G')
                else:
                    sol[i], sol[mi] = ('U', 'A') if bp_dir == 0 else ('A', 'U')

        # Long-range pairs: AUC-limited alphabet
        if target[i] == '(' and ci == 'N' and cm == 'N' and (mi - i) > LONGRANGE:
            bp_type = random.randint(0, 2)
            pairs_lr = [('A', 'U'), ('C', 'G'), ('U', 'A')]
            sol[i], sol[mi] = pairs_lr[bp_type]

        # '{' dovetail hints: force 100% GC
        if init_target[i] == '{' and ci == 'N':
            if cm == 'A':
                sol[i] = 'U'; sol[mi] = 'A'
            elif cm == 'U':
                sol[i] = 'A'; sol[mi] = 'U'
            elif cm == 'C':
                sol[i] = 'G'; sol[mi] = 'C'
            elif cm == 'G':
                sol[i] = 'C'; sol[mi] = 'G'
            elif cm == 'N':
                bp_dir = random.randint(0, 1)
                sol[i], sol[mi] = ('G', 'C') if bp_dir == 0 else ('C', 'G')

    return sol

def append_trace_analysis(output_path: str, seq: str) -> None:
    """
    Match dragon.pl preview():
    - write seq.txt
    - run trace_analysis.py pattern.txt seq.txt
    - append stdout to the exported design file
    """
    clean_seq = ''.join(ch for ch in seq.strip().upper() if ch in "AUCGN")

    with open("seq.txt", "w", encoding="utf-8") as f:
        f.write(clean_seq)

    with open(output_path, "a", encoding="utf-8") as fout:
        subprocess.run(
            [sys.executable, "trace_analysis.py", "pattern.txt", "seq.txt"],
            stdout=fout,
            stderr=subprocess.STDOUT,
            check=False,
        )

# ============================================================
# Export / Output
# ============================================================

def export_result(best_sol: List[str], name: str, output_file_path: str,
                  target_pk_str: str, gc_stats: dict, kl_analysis: dict,
                  failed_initial_design: bool, failed_pattern_test: bool,
                  lowest_ed: float, log) -> tuple:
    """
    Compute ensemble diversity; export design file if it is a new optimum.
    Returns (lowest_ed, keep_ed).
    Fixes the broken qr// regex in the original sub export.
    """
    seq = ''.join(best_sol)
    ed = rna_ensemble_diversity(seq)
    log(f"Ensemble Diversity = {ed:.4f}\n\n")

    keep_ed = False
    if ed < lowest_ed:
        lowest_ed = ed

        if failed_initial_design:
            log(" WARNING: Misfolding within the locked-sequence required altering the target structure.\n")
            log("          Sequence design failed for the inputted target structure.\n\n")

        log(f"NEW OPTIMAL: Exporting design to {name}_design\n\n{seq}\n\n")

        filename_tag = int(lowest_ed * 100) / 100
        out_path = f"{output_file_path}{name}_design_{filename_tag}.txt"
        with open(out_path, 'w', encoding='utf-8') as f:
            f.write(f"{name}\n")
            f.write(f"{target_pk_str}\n")
            f.write(f"{seq}\n\n")

            f.write(
                f" GC Content: {gc_stats['gc_cont']}  TotalPairs: {gc_stats['pairs']}\n"
                f" GC pairs: {gc_stats['gc_pairs']}     AU pairs: {gc_stats['au_pairs']}"
                f"     GU pairs: {gc_stats['gu_pairs']}"
                f"    KL: {gc_stats['kl_gc_cont']} %GC\n\n\n"
            )

            f.write(f"Ensemble Diversity: {lowest_ed} \n\n")

            if failed_initial_design:
                f.write(" WARNING: Misfolding within the locked-sequence required "
                        "altering the target structure.\n")
                f.write("          Sequence design failed for the inputted target structure.\n\n")

            if failed_pattern_test:
                f.write(" WARNING: Constrained sequences contain undesired seqeunce "
                        "patternds/duplication. \n")
                f.write("          Check pattern map carefully before using this design.\n\n")

            f.write("Kissing Loop List:\n")
            kl = kl_analysis
            for j in range(1, kl.get('num_kl', 0) + 1):
                f.write(
                    f"{j} A,{kl['a_list'][j]},"
                    f"{j} B,{kl['b_list'][j]},"
                    f"{kl['e_list'][j]}\n"
                )

            f.write(f"\n{kl.get('non_cognates', '')}")
            f.write("\n\n\n")

        append_trace_analysis(out_path, seq)
        keep_ed = True
        
    return lowest_ed, keep_ed

# ============================================================
# Banners
# ============================================================

def print_welcome(log):
    log("\n____    ____ _  _ ____ _    _  _ ____  \n")
    log("|__/ __ |___ |  | |  | |    |  | |__/ \n")
    log("|  \\    |___  \\/  |__| |___  \\/  |  \\  v2.03\n\n\n")
    log(" Dragon Version - Ensemble Diversity Optimizer (Python Port)\n")
    log("                                                  \n")
    log("       0                                                                                               \n")
    log("      1111111                             1111  0                 1    0                            111\n")
    log("     11 010111                        11111000 0              11111 1111                      1    011\n")
    log("  0111  1111110110                 111111111   0110          00 0 111111110                   000111 01\n")
    log("1 10111 1110  11111              000011111111 101001 1    00  1000 11111100111110            10 011100111\n")
    log("11011001101 111111110           00001     111 100 01000  111100 00111     0111110       000011   1  11\n")
    log("10101 1000  1111100110         0000 00    1 11111 11101  011100001 11     0 111111    010101      11   \n")
    log("1000 1110   1110 011111      0000   00     1 1111111111  111 000101 1     00  1110111 0110000      11   \n")
    log("010 0000    1     01111     0000   000     11 111111000 1110 0010  1     000    1011110010000 1     11  \n")
    log("00000001   11      10111   000      00     1  111111000 1110000  11     00       0111111    000     1  \n")
    log("010011     1        10111 00001      00     1   0111000 011111   1     000       001111      00        \n")
    log(" 10001    11    1111 101111111111      0         01111  11111           0       01010111     00     00 \n")
    log("1000      11     111 1 0110111 111             0  0111011011    0             01010110111    000     0 \n")
    log("001       1     1111 0 11010 11 11     1     000   0111000     00             0000011100101  00      0 \n")
    log("011                1  0001000 0        11     00  1101110 11  000     1     101111010111000  00      00\n")
    log("11        0     0    0000 0010   11     11    000 1000111010 000     11     00111 111  11000          0\n")
    log("        00      0   00100 10100   1     11      000000001111 00      1     11001  0111 011011   11    0\n")
    log("       000     00 000000110001110 11     1       00001 111100       11    1111111 0111 0 01110  11     \n")
    log("1      00      0 00000111011001110 1     11     00000   001010      1    11111 111 1110111111000  1     \n")
    log("110   0        1001  100001111 011111     111110000      000010    11   1111  111  1100001110001  1   1\n")
    log(" 11101      11100 1  11 1 1  1   111 1     1100000         010111111  11100   111       11  1111  11111\n")
    log("  111110011110 0  1 1 1   0  00  111111111100001             11110011110001  111        11   0111111111\n")
    log("      0100001   11           100 1111111100                    101010  101111            1  001 111011\n")
    log("         11111   1               0  111111                         101011 0                 1101110    \n")
    log("         0                       0     1                             010 0                   0 110     \n")
    log("\n\nWritten by Cody Geary. Copyright 2023. Python port 2025.\n")


def print_victory(log, start_time: float, failed_initial_design: bool):
    log("                                 .''.                    \n")
    log("       .''.             *''*    :_\\/_:     .            \n")
    log("      :_\\/_:   .    .:.*_\\/_*   : /\\ :  .'.:.'.       \n")
    log("  .''.: /\\ : _\\(/_  ':'* /\\ *  : '..'.  -=:o:=-       \n")
    log(" :_\\/_:'.:::. /)\\*''*  .|.* '.\\'/.'_\\(/_'.':'.'      \n")
    log(" : /\\ : :::::  '*_\\/_* | |  -= o =- /)\\    '  *       \n")
    log("  '..'  ':::'   * /\\ * |'|  .'/.\\'.  '._____           \n")
    log("      *        __*..* |  |     :      |.   |' .---'|     \n")
    log("       _*   .-'   '-. |  |     .--'|  ||   | _|    |     \n")
    log("    .-'|  _.|  |    ||   '-__  |   |  |    ||      |     \n")
    log("    |' | |.    |    ||       | |   |  |    ||      |     \n")
    log(" ___|  '-'     '    ''       '-'   '-.'    '`      |____ \n")
    log("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n\n\n")
    elapsed = time.time() - start_time
    log(f"Total Design Time: {elapsed:.1f} seconds\n\n")
    if failed_initial_design:
        log(" WARNING: Misfolding within the locked-sequence required altering the target structure.\n")
        log("          Sequence design failed for the inputted target structure.\n\n")

# ============================================================
# Main
# ============================================================

def main():
    start_time = time.time()

    # --- Parse CLI arguments ---
    output_file_path = sys.argv[1] if len(sys.argv) > 1 else None

    # --- Read target file ---
    name, target_pk_str, constraint_str, start_seq_str = read_target('target.txt')

    # --- Parse structure into three representations ---
    init_target, target_nopk, target_scrubbed = parse_target(target_pk_str)
    target = target_scrubbed  # main working copy (normalised brackets)
    target_str = ''.join(target_nopk)  # for bp_distance (no pseudoknots)
    strand_length = len(target_str)

    # --- Build pair map from original PK notation ---
    pair_map = build_pair_map(target_pk_str)

    # --- Validate constraint ---
    constraint = list(constraint_str)
    if len(constraint) != strand_length:
        print("\nInvalid Constraint, using all N instead.")
        constraint = ['N'] * strand_length

    # --- Setup output directory ---
    if output_file_path is None:
        os.makedirs(name, exist_ok=True)
        output_file_path = name
    output_file_path = output_file_path.rstrip('/') + '/'

    spool_path = f"{output_file_path}{name}_spool.txt"
    spool_file = open(spool_path, 'w')

    def log(msg: str):
        """Write to spool and stdout simultaneously."""
        spool_file.write(msg)
        print(msg, end='', flush=True)

    def log_gc(sol):
        gs = count_gc(sol, target, init_target, pair_map, strand_length)
        log(f"\n GC Content: {gs['gc_cont']}    TotalPairs: {gs['pairs']}\n"
            f" GC pairs: {gs['gc_pairs']}     AU pairs: {gs['au_pairs']}"
            f"     GU pairs: {gs['gu_pairs']}    KL: {gs['kl_gc_cont']} %GC\n")
        return gs

    log(f"Output path is {output_file_path}\n")
    print_welcome(log)

    # --- Generate or load initial sequence ---
    generate_sequence = True
    if start_seq_str:
        log("I found a sequence! Inserting it.\n")
        best_sol = list(start_seq_str)
        generate_sequence = False
    else:
        best_sol = initialize_sequence(target, init_target, constraint,
                                       pair_map, strand_length)

    log(f"Designing {name}:\n{target_pk_str}\n")
    log(f"The design is {strand_length} nts long\n\n")
    log(f"Target:\n{target_str}\n\nConstraint\n{constraint_str}\n\n\n")

    # --- Fold memoisation cache ---
    fold_cache: dict = {}

    def cached_fold(seq: str) -> tuple:
        """Return (structure, previously_tried_bool)."""
        if seq in fold_cache:
            log("<found cached fold>\n")
            return fold_cache[seq], True
        struct = rna_fold(seq)
        fold_cache[seq] = struct
        return struct, False

    # --- Initial fold ---
    best_sol_seq = ''.join(best_sol)
    parent_fold_str, _ = cached_fold(best_sol_seq)
    parent_fold = list(parent_fold_str)
    parent_dist = rna_distance(parent_fold_str, target_str)

    log(f"Initial Sequence: Distance = {parent_dist}\n")
    gc_stats = log_gc(best_sol)
    log(f"{best_sol_seq}\n{parent_fold_str}\n")

    trial_sol = list(best_sol)
    mut_map   = ['-'] * strand_length
    defect_map = ['-'] * strand_length
    rad_level  = FAV_RAD_LEVEL
    drift_rate = MY_DRIFT_RATE
    drift_period = DRIFT_PERIOD
    attempts   = 0
    distance   = parent_dist
    failed_initial_design = False
    failed_pattern_test   = False
    init_setup = INIT_SETUP

    # ============================================================
    # PHASE 1: Initial setup loop
    # Coarsely assigns pairs/loops to approach the target structure
    # ============================================================

    while attempts < init_setup:
        if parent_dist == 0 or not generate_sequence:
            break
        attempts += 1
        helix_end = 0

        for i in range(1, strand_length - 1):
            ci  = constraint[i]
            mi  = pair_map[i]
            cmi = constraint[mi] if mi != i else 'N'

            if ci == 'N' and cmi == 'N':
                bp_dir  = random.randint(0, 1)
                bp_type = random.uniform(0, 100)

                if target[i] == '.' and parent_fold[i] != '.':
                    trial_sol[i] = 'A'

                if target[i] in ('(', '[') and i + 1 < strand_length and target[i+1] not in ('(', '['):
                    helix_end = 1
                if target[i] in (')', ']') and i + 1 < strand_length and target[i+1] not in (')', ']'):
                    helix_end = 1
                if helix_end:
                    bp_type = 100

                if target[i] == '(' and parent_fold[i] != '(':
                    if bp_type < (100 - MAX_GC - TARGET_GU):
                        trial_sol[i], trial_sol[mi] = ('A', 'U') if bp_dir == 0 else ('U', 'A')
                    else:
                        trial_sol[i], trial_sol[mi] = ('G', 'C') if bp_dir == 0 else ('C', 'G')

                if (i + 1 < strand_length and target[i] == ')' and target[i+1] == '(' and
                        parent_fold[i] != ')' and parent_fold[i+1] != '('):
                    trial_sol[mi] = 'G'; trial_sol[i] = 'C'
                    trial_sol[pair_map[i+1]] = 'C'; trial_sol[i+1] = 'G'

                if (i > 0 and target[i] == '(' and target[i-1] == ')' and
                        parent_fold[i] != '(' and parent_fold[i-1] != ')'):
                    trial_sol[mi] = 'C'; trial_sol[i] = 'G'
                    trial_sol[pair_map[i-1]] = 'G'; trial_sol[i-1] = 'C'

                if (i + 1 < strand_length and target[i] == '(' and target[i+1] == '(' and
                        (mi - pair_map[i+1]) > 1 and
                        parent_fold[i] != '(' and parent_fold[i+1] != '('):
                    trial_sol[mi] = 'G'; trial_sol[i] = 'C'
                    trial_sol[pair_map[i+1]] = 'G'; trial_sol[i+1] = 'C'

                if (i + 1 < strand_length and target[i] == ')' and target[i+1] == ')' and
                        (mi - pair_map[i+1]) > 1 and
                        parent_fold[i] != ')' and parent_fold[i+1] != ')'):
                    trial_sol[mi] = 'G'; trial_sol[i] = 'C'
                    trial_sol[pair_map[i+1]] = 'G'; trial_sol[i+1] = 'C'

            elif ci == 'S' and cmi == 'S':
                if (i + 1 < strand_length and target[i] == ')' and target[i+1] == '(' and
                        parent_fold[i] != ')' and parent_fold[i+1] != '('):
                    trial_sol[mi] = 'G'; trial_sol[i] = 'C'
                    trial_sol[pair_map[i+1]] = 'C'; trial_sol[i+1] = 'G'
                if (i > 0 and target[i] == '(' and target[i-1] == ')' and
                        parent_fold[i] != '(' and parent_fold[i-1] != ')'):
                    trial_sol[mi] = 'C'; trial_sol[i] = 'G'
                    trial_sol[pair_map[i-1]] = 'G'; trial_sol[i-1] = 'C'
                if (i + 1 < strand_length and target[i] == '(' and target[i+1] == '(' and
                        (mi - pair_map[i+1]) > 1 and
                        parent_fold[i] != '(' and parent_fold[i+1] != '('):
                    trial_sol[mi] = 'G'; trial_sol[i] = 'C'
                    trial_sol[pair_map[i+1]] = 'G'; trial_sol[i+1] = 'C'
                if (i + 1 < strand_length and target[i] == ')' and target[i+1] == ')' and
                        (mi - pair_map[i+1]) > 1 and
                        parent_fold[i] != ')' and parent_fold[i+1] != ')'):
                    trial_sol[mi] = 'G'; trial_sol[i] = 'C'
                    trial_sol[pair_map[i+1]] = 'G'; trial_sol[i+1] = 'C'

        new_trial_seq = ''.join(trial_sol)
        new_fold_str, _ = cached_fold(new_trial_seq)
        trial_dist = rna_distance(new_fold_str, target_str)

        if trial_dist < parent_dist * MISFOLD_TOLERANCE:
            for i in range(strand_length):
                mut_map[i]    = trial_sol[i] if best_sol[i] != trial_sol[i] else '-'
                defect_map[i] = 'X' if parent_fold[i] != target_nopk[i] else '_'
                best_sol[i]   = trial_sol[i]

            best_sol_seq    = ''.join(trial_sol)
            parent_dist     = trial_dist
            parent_fold_str = new_fold_str
            parent_fold     = list(parent_fold_str)
            distance        = trial_dist

            log(f"\nIteration# {attempts}, Distance = {trial_dist}")
            gc_stats = log_gc(best_sol)
            if attempts > 1:
                log(''.join(mut_map) + '\n')
            log(f"{best_sol_seq}\n{parent_fold_str}\n{''.join(defect_map)}\n")
        else:
            log(f"( {attempts} )")
            trial_sol = list(best_sol)

    # ============================================================
    # PHASE 2: Directed mutation cycle
    # Targeted mutations fix structural mismatches; periodic drift
    # explores sequence space.  Bug fix: drift_event_count now only
    # resets when a drift event actually fires (not every iteration).
    # ============================================================

    log("\n__________________________________________________________\n")
    log("Begin Directed Mutation Cycle\n")

    drift_event_count    = 0
    improved_parent_count = 0
    constraint_check     = 0
    failed_tries         = 0

    while distance > 0:
        attempts            += 1
        improved_parent_count += 1
        drift_event_count   += 1

        mutation_sites     = [0] * strand_length
        num_mutation_spots = 0
        is_drift           = (drift_event_count == drift_period)

        log(f"({attempts})")
        log("drift->" if is_drift else "target->")
        if improved_parent_count > 1 and not is_drift:
            log("expanded->")

        for i in range(2, strand_length - 2):
            ci  = constraint[i]
            mi  = pair_map[i]

            mutate_now = False
            if is_drift and random.uniform(0, 100) < drift_rate and ci in ('N', 'K', 'S'):
                mutate_now = True

            if (parent_fold[i] != target_nopk[i] and not is_drift and
                    random.uniform(0, 100) < SITE_MUTATION_RATE and ci in ('N', 'K', 'S')):
                mutate_now = True

            if improved_parent_count > 1:
                if (not is_drift and ci == 'N' and
                        random.uniform(0, 100) < SITE_MUTATION_RATE and
                        (parent_fold[i+1] != target_nopk[i+1] or
                         parent_fold[i-1] != target_nopk[i-1] or
                         parent_fold[i+2] != target_nopk[i+2] or
                         parent_fold[i-2] != target_nopk[i-2])):
                    mutate_now = True

                if parent_fold[i] != target[i] and target[i] == ')' and target[i+1] == '(':
                    for idx in [i, mi, i+1, pair_map[i+1]]:
                        if 0 <= idx < strand_length:
                            mutation_sites[idx] = 1
                    if mi > 0 and mi - 1 >= 0:
                        mutation_sites[mi - 1] = 1
                        mm = pair_map[mi - 1]
                        if 0 <= mm < strand_length:
                            mutation_sites[mm] = 1
                    mutate_now = True

                for offset in [0, 1, -1]:
                    idx = pair_map[i + offset] if 0 <= i + offset < strand_length else -1
                    if 0 <= idx < strand_length and parent_fold[idx] != target[idx]:
                        mutation_sites[idx] = 1
                        mutate_now = True

            if mutate_now:
                mutation_sites[i] = 1
                num_mutation_spots += 1

        mut_scale = rad_level / (num_mutation_spots + 0.001)
        log(f" M[{num_mutation_spots}] (*{rad_level}) ")

        for i in range(2, strand_length - 2):
            if random.uniform(0, 1) >= mut_scale or mutation_sites[i] != 1:
                continue
            ci  = constraint[i]
            mi  = pair_map[i]
            cmi = constraint[mi] if mi != i else 'N'
            ti  = target[i]

            if ti in ('(', ')') and ci == 'N' and cmi == 'N':
                log(f" {i}")
                coin = random.randint(0, 1)

                if random.uniform(0, 100) < DOUBLE_MUTANT:
                    log(":")
                    if random.uniform(0, 100) < GC_MUT_RATIO:
                        if coin == 1:
                            trial_sol[i] = 'G'
                            trial_sol[mi] = 'C'
                        else:
                            trial_sol[i] = 'C'
                            trial_sol[mi] = 'G'
                    else:
                        if coin == 1:
                            trial_sol[i] = 'A'
                            trial_sol[mi] = 'U'
                        else:
                            trial_sol[i] = 'U'
                            trial_sol[mi] = 'A'
                else:
                    log(".")
                    if best_sol[i] == 'G' and best_sol[mi] == 'C':
                        trial_sol[mi] = 'U'
                    elif best_sol[i] == 'C' and best_sol[mi] == 'G':
                        trial_sol[i] = 'U'
                    elif best_sol[i] == 'A' and best_sol[mi] == 'U':
                        trial_sol[i] = 'G'
                    elif best_sol[i] == 'U' and best_sol[mi] == 'A':
                        trial_sol[mi] = 'G'
                    elif best_sol[i] == 'G' and best_sol[mi] == 'U':
                        if coin == 0:
                            trial_sol[mi] = 'C'
                        else:
                            trial_sol[i] = 'A'
                    elif best_sol[i] == 'U' and best_sol[mi] == 'G':
                        if coin == 0:
                            trial_sol[i] = 'C'
                        else:
                            trial_sol[mi] = 'A'
            elif ti == '.' and ci == 'N':
                log("L")
                if random.uniform(0, 100) < MUTATE_LOOP:
                    trial_sol[i] = rand_base()
            elif ti == '.' and ci == 'R':
                log("Lr"); trial_sol[i] = rand_purine()
            elif ti == '.' and ci == 'Y':
                log("Lr"); trial_sol[i] = rand_pyrimidine()
            elif ti == '.' and ci == 'S':
                log("Lr"); trial_sol[i] = rand_strong()
            elif ti == '.' and ci == 'W':
                log("Lr"); trial_sol[i] = rand_weak()

            if ti in ('[', ']') and ci == 'N' and cmi == 'N':
                log(f" {i}[:")
                coin = random.randint(0, 1)
                if random.uniform(0, 100) < KL_GC:
                    trial_sol[i], trial_sol[mi] = ('G', 'C') if coin else ('C', 'G')
                else:
                    trial_sol[i], trial_sol[mi] = ('A', 'U') if coin else ('U', 'A')

            if ti in ('(', '[', ')', ']') and ci == 'K' and cmi == 'K':
                log("K")
                if best_sol[i] == 'G':
                    trial_sol[mi] = 'G'; trial_sol[i] = 'U'
                elif best_sol[i] == 'U':
                    trial_sol[i] = 'G'; trial_sol[mi] = 'U'

            if ti in ('(', '[', ')', ']') and ci == 'S' and cmi == 'S':
                log("S")
                if best_sol[i] == 'G':
                    trial_sol[mi] = 'G'; trial_sol[i] = 'C'
                elif best_sol[i] == 'C':
                    trial_sol[i] = 'G'; trial_sol[mi] = 'C'

        new_trial_seq = ''.join(trial_sol)
        new_fold_str, previously_tried = cached_fold(new_trial_seq)
        new_fold_array = list(new_fold_str)
        trial_dist = rna_distance(new_fold_str, target_str)
        distance   = trial_dist

        if trial_dist < parent_dist:
            improved_parent_count = 0
        keep_sequence = trial_dist < parent_dist * MISFOLD_TOLERANCE

        if keep_sequence and not previously_tried:
            rad_level += 1
            for i in range(strand_length):
                mut_map[i]    = trial_sol[i] if best_sol[i] != trial_sol[i] else '-'
                defect_map[i] = 'X' if new_fold_array[i] != target_nopk[i] else '-'
                best_sol[i]   = trial_sol[i]

            best_sol_seq = ''.join(best_sol)
            constraint_check = sum(
                1 for i in range(strand_length)
                if defect_map[i] == 'X' and constraint[i] in ('N', 'K', 'S')
            )
            display = trial_dist < parent_dist or OUTPUT_THRESHOLD == 1 or constraint_check == 0

            parent_dist     = trial_dist
            parent_fold_str = new_fold_str
            parent_fold     = list(parent_fold_str)

            if display:
                log(f"\n\nIteration# {attempts}, Distance = {trial_dist}\n")
                gc_stats = log_gc(best_sol)
                log(f"\nMutations:\n{''.join(mut_map)}\n{best_sol_seq}\n{parent_fold_str}\n")
                log(f"Defect Map:\n{''.join(defect_map)}\n")
            else:
                log(f" keep -D{trial_dist}")

            drift_event_count = 0

        else:
            rad_level = max(1, rad_level - 1)
            log(f" dump -D{trial_dist}")
            trial_sol = list(best_sol)

        if previously_tried:
            drift_event_count = drift_period - 1
            log("Forcing Drift")
            rad_level += 1

        log("\n")

        if constraint_check == 0 and distance > 0 and keep_sequence:
            failed_tries += 1
            drift_period = 1
            log(f" <Blocked by Constraints - {failed_tries} >\n")

        if (constraint_check == 0 and distance > 0 and
                failed_tries > MAX_FAILED_TRIES and keep_sequence):
            log("All defects are locked by sequence constraints. "
                "Resetting target to current best defect mask\n")
            distance = 0
            for i in range(strand_length):
                if target[i] not in ('[', ']'):
                    target[i] = new_fold_array[i]
                target_nopk[i] = new_fold_array[i]
            target_str = new_fold_str.strip()
            failed_initial_design = True
            log(f"WARNING: Target reset to best-effort structure after {MAX_FAILED_TRIES} attempts.\n")

    log(f"Success!\nDesign {name} folded.\n")

    _, kl_repeats_val = pattern_prevent(trial_sol, target_nopk, constraint, strand_length)
    log(f"\n{kl_repeats_val} KL repeats found.\n\n")

    # ============================================================
    # PHASE 3: GC-Reduction & Pattern Cleanup
    # ============================================================

    log("\n__________________________________________________________\n")
    log("Begin GC-Reduction & GC/UA-Rich Reduction\n")

    rad_level  = FAV_RAD_LEVEL
    drift_rate = MY_DRIFT_RATE

    for i in range(strand_length):
        trial_sol[i] = best_sol[i]

    rpt = count_repeats(trial_sol, target, constraint, strand_length,
                        COMPLEMENT_WINDOW, DUPLICATE_WINDOW, PALINDROME_WINDOW)
    last_pr  = rpt['pattern_repeats']
    last_poly = rpt['poly_repeats']
    last_rs  = rpt['restriction_sites']
    last_cz  = rpt['complement_zones']
    last_dz  = rpt['duplication_zones']
    last_pal = rpt['palindromes']

    log(f"\nThere are {last_cz} Complementary Regions, {last_pr} Strong/Weak Regions, "
        f"{last_poly} Poly-N5 Repeats, {last_dz} nt Duplicated sequence, "
        f"and {last_rs} common restriction sites present to start:\n")

    def _pattern_score(pr, poly, rs, cz, dz, pal):
        return pr + poly + rs + cz / COMPLEMENT_WINDOW + dz + pal

    def _all_zero(pr, poly, rs, cz, dz, pal):
        return pr == 0 and poly == 0 and rs == 0 and cz == 0 and dz == 0 and pal == 0

    def _is_masked(rpt_map):
        """Return (unmasked_count, masked_count) of pattern positions."""
        unmasked = sum(1 for idx in range(strand_length)
                       if constraint[idx] in ('N', 'S', 'K', 'R', 'Y') and rpt_map[idx] != '-')
        masked   = sum(1 for idx in range(strand_length)
                       if constraint[idx] in ('A', 'U', 'C', 'G') and rpt_map[idx] != '-')
        return unmasked, masked

    happy      = False
    test_first = True
    gu_put     = [0] * strand_length

    while not happy:
        rpt = count_repeats(trial_sol, target, constraint, strand_length,
                            COMPLEMENT_WINDOW, DUPLICATE_WINDOW, PALINDROME_WINDOW)
        repeat_map    = rpt['repeat_map']
        pr_cnt        = rpt['pattern_repeats']
        poly_cnt      = rpt['poly_repeats']
        rs_cnt        = rpt['restriction_sites']
        cz_cnt        = rpt['complement_zones']
        dz_cnt        = rpt['duplication_zones']
        pal_cnt       = rpt['palindromes']

        if not test_first:
            attempts += 1

            # Inner retry: keep trying until patterns don't worsen
            while True:
                log(f"(Iteration {attempts})")
                num_mutation_spots = 0
                mutation_sites = [0] * strand_length

                for i in range(strand_length):
                    trial_sol[i] = best_sol[i]
                    gu_put[i]    = 0
                    mutation_sites[i] = 0

                gs = count_gc(best_sol, target, init_target, pair_map, strand_length)

                for i in range(strand_length):
                    mutate_now = False
                    coin = random.randint(0, 1)
                    rand1 = random.uniform(0, 100)

                    if 1 < i < strand_length - 1:
                        # Random GU placement
                        if (random.randint(0, max(1, int(100 / (TARGET_GU / 10 + 1)))) == 0 and
                                gs['gu_pairs'] < TARGET_GU * gs['pairs'] * 0.01 and
                                target[i] == '('):
                            gu_put[i] = 1; mutate_now = True
                            mutation_sites[i] += 1

                        if (repeat_map[i] == 'P' and random.uniform(0, 100) < 15 and
                                not mutate_now and target[i] == '('):
                            gu_put[i] = 1; mutate_now = True
                            mutation_sites[i] += 1

                        # Suppress GU adjacent to another GU
                        if target[i] in ('(', ')') and mutate_now:
                            adj_gu = False
                            for nb in [i - 1, i + 1]:
                                if 0 <= nb < strand_length and pair_map[nb] != nb:
                                    a2, b2 = trial_sol[nb], trial_sol[pair_map[nb]]
                                    if {a2, b2} == {'G', 'U'}:
                                        adj_gu = True
                            if adj_gu:
                                gu_put[i] = 0; mutate_now = False
                                mutation_sites[i] = max(0, mutation_sites[i] - 1)

                    rm = repeat_map[i]
                    if rm == '-':
                        if target[i] in ('(', '[') and rand1 < drift_rate / 2:
                            mutate_now = True; mutation_sites[i] += 1
                    elif rm == 'U':
                        mutate_now = True; mutation_sites[i] += 2
                    elif rm in ('G', 'C', 'A'):
                        mutate_now = True; mutation_sites[i] += 1
                    elif rm in ('S', 'W', 'X', 'D'):
                        if coin == 0:
                            mutate_now = True; mutation_sites[i] += 1
                    elif rm == 'P':
                        if random.uniform(0, max(1, cz_cnt)) < 5:
                            mutate_now = True; mutation_sites[i] += 1
                    elif rm == 'L':
                        if random.uniform(0, max(1, pal_cnt)) < 5:
                            mutate_now = True; mutation_sites[i] += 1

                    if mutate_now:
                        num_mutation_spots += 1

                rad_level  = max(2, min(20, rad_level))
                mut_scale  = rad_level / (num_mutation_spots + 0.001)
                log(f" M[{num_mutation_spots}] (*{rad_level}) ")

                for i in range(strand_length - 1):
                    if random.uniform(0, 1) >= mut_scale * mutation_sites[i]:
                        continue
                    ci  = constraint[i]
                    mi  = pair_map[i]
                    cmi = constraint[mi] if mi != i else 'N'
                    ti  = target[i]
                    rm  = repeat_map[i]
                    coin = random.randint(0, 1)

                    if ti == '.' and ci == 'N':
                        log("."); trial_sol[i] = rand_base()

                    if rm == 'P' and gu_put[i] == 1:
                        if ci == 'N' and cmi == 'N':
                            log("+W ")
                            trial_sol[i], trial_sol[mi] = ('G', 'U') if coin == 0 else ('U', 'G')
                        elif ci == 'R' and cmi in ('Y', 'N'):
                            log("+W "); trial_sol[i] = 'G'; trial_sol[mi] = 'U'
                        elif ci == 'Y' and cmi in ('R', 'N'):
                            log("+W "); trial_sol[i] = 'U'; trial_sol[mi] = 'G'

                    if (rm in ('G', 'C') and random.uniform(0, 1) < 0.5 and
                            ti in ('(', ')', '[', ']') and ci == 'N' and cmi == 'N'):
                        log(":fS ")
                        trial_sol[i], trial_sol[mi] = ('G', 'C') if trial_sol[i] == 'C' else ('C', 'G')

                    elif (ti in ('(', ')', '[', ']') and
                          ((ci == 'R' and cmi in ('Y', 'N')) or
                           (ci == 'Y' and cmi in ('R', 'N')))):
                        if (rm in ('G', 'A') and ci == 'R' and
                                random.uniform(0, 1) < 0.5 and cmi in ('N', 'Y')):
                            log(":R ")
                            trial_sol[i], trial_sol[mi] = ('A', 'U') if trial_sol[i] == 'G' else ('G', 'C')

                        if (rm in ('X', 'U', 'C') and ci == 'Y' and
                                random.uniform(0, 1) < 0.5 and cmi in ('N', 'R')):
                            log(":Y ")
                            trial_sol[i], trial_sol[mi] = ('C', 'G') if trial_sol[i] == 'U' else ('U', 'A')

                    elif (ti in ('(', ')', '[', ']') and ci == 'N' and cmi == 'N'):
                        log(":")

                        # GC -> AU
                        if ((trial_sol[i] == 'C' and trial_sol[mi] == 'G') or
                            (trial_sol[i] == 'G' and trial_sol[mi] == 'C')):
                            gc_pct = gs['gc_pairs'] / max(1, gs['pairs']) * 100
                            if gc_pct > MAX_GC or random.uniform(0, 1) < 0.5:
                                log("+A ")
                                if coin == 0:
                                    trial_sol[i], trial_sol[mi] = 'A', 'U'
                                else:
                                    trial_sol[i], trial_sol[mi] = 'U', 'A'

                        # AU -> GC
                        elif ((trial_sol[i] == 'A' and trial_sol[mi] == 'U') or
                              (trial_sol[i] == 'U' and trial_sol[mi] == 'A')):
                            gc_pct = gs['gc_pairs'] / max(1, gs['pairs']) * 100
                            if gc_pct < MAX_GC or random.uniform(0, 1) < 0.25:
                                log("+G ")
                                if coin == 0:
                                    trial_sol[i], trial_sol[mi] = 'G', 'C'
                                else:
                                    trial_sol[i], trial_sol[mi] = 'C', 'G'

                        # GU flip / delete
                        elif ((trial_sol[i] == 'G' and trial_sol[mi] == 'U') or
                              (trial_sol[i] == 'U' and trial_sol[mi] == 'G')):
                            if random.uniform(0, 1) < 0.9:
                                log(":fW ")
                                if trial_sol[i] == 'G':
                                    trial_sol[i], trial_sol[mi] = 'U', 'G'
                                else:
                                    trial_sol[i], trial_sol[mi] = 'G', 'U'
                            else:
                                log(":xW ")
                                if trial_sol[i] == 'G':
                                    trial_sol[i] = 'A'
                                else:
                                    trial_sol[mi] = 'A'

                    elif (ti in ('(', ')', '[', ']') and ci == 'K' and cmi == 'K'):
                        log(":fW ")
                        if trial_sol[i] == 'G':
                            trial_sol[i], trial_sol[mi] = 'U', 'G'
                        else:
                            trial_sol[i], trial_sol[mi] = 'G', 'U'

                    elif (ti in ('(', ')', '[', ']') and ci == 'S' and cmi == 'S'):
                        log(":fS ")
                        if trial_sol[i] == 'G':
                            trial_sol[i], trial_sol[mi] = 'C', 'G'
                        else:
                            trial_sol[i], trial_sol[mi] = 'G', 'C'

                rpt2 = count_repeats(trial_sol, target, constraint, strand_length,
                                     COMPLEMENT_WINDOW, DUPLICATE_WINDOW, PALINDROME_WINDOW)
                pr_cnt   = rpt2['pattern_repeats']
                poly_cnt = rpt2['poly_repeats']
                rs_cnt   = rpt2['restriction_sites']
                cz_cnt   = rpt2['complement_zones']
                dz_cnt   = rpt2['duplication_zones']
                pal_cnt  = rpt2['palindromes']
                repeat_map = rpt2['repeat_map']

                log(f"(SW{pr_cnt}/N{poly_cnt}/R{rs_cnt}/W{cz_cnt}/L{pal_cnt})")

                if (_pattern_score(pr_cnt, poly_cnt, rs_cnt, cz_cnt, dz_cnt, pal_cnt) >
                        _pattern_score(last_pr, last_poly, last_rs, last_cz, last_dz, last_pal)):
                    log(" - bad roll, trying again...\n")
                    rad_level -= 1
                else:
                    break  # Accept this mutation round

        test_first = False

        unmasked, masked = _is_masked(repeat_map)
        if unmasked == 0 and masked > 0:
            log("\n\n Remaining patterns are masked by sequence constraints.  Passing to next phase. \n\n")
            happy = True

        new_trial_seq = ''.join(trial_sol)
        log("SUCCESS folding...")
        new_fold_str, _ = cached_fold(new_trial_seq)
        trial_dist = rna_distance(new_fold_str, target_str)

        keep_sequence = (trial_dist == 0 and
                         (pr_cnt <= last_pr or poly_cnt <= last_poly or
                          rs_cnt <= last_rs or cz_cnt > last_cz or dz_cnt > last_dz))

        if keep_sequence:
            last_pr = pr_cnt; last_poly = poly_cnt; last_rs = rs_cnt
            last_cz = cz_cnt; last_dz = dz_cnt; last_pal = pal_cnt
            rad_level += 2
            for i in range(strand_length):
                mut_map[i]  = trial_sol[i] if best_sol[i] != trial_sol[i] else '-'
                best_sol[i] = trial_sol[i]

            best_sol_seq    = ''.join(best_sol)
            parent_dist     = trial_dist
            parent_fold_str = new_fold_str
            parent_fold     = list(parent_fold_str)
            gc_stats        = log_gc(best_sol)

            log(f" 8S|8W: {pr_cnt}  5N: {poly_cnt}  Restrict: {rs_cnt}  "
                f"Complement: {cz_cnt}  Dup: {dz_cnt}  Palindromes: {pal_cnt}\n")
            log(f"Mutation Map:\n{''.join(mut_map)}\n{best_sol_seq}\n")
            log(f"Patterns Remaining:\n{''.join(repeat_map)}\n\n")

            unmasked, masked = _is_masked(repeat_map)
            if unmasked == 0 and masked > 0:
                log("\n\n Remaining patterns are masked by sequence constraints.  Passing to next phase. \n\n")
                happy = True
        else:
            rad_level = max(3, rad_level - 1)
            log(f" Bad Fold: -D{trial_dist}")
            trial_sol = list(best_sol)

        gc_stats = count_gc(best_sol, target, init_target, pair_map, strand_length)
        if (gc_stats['gc_cont'] < MAX_GC and _all_zero(pr_cnt, poly_cnt, rs_cnt, cz_cnt, dz_cnt, pal_cnt)
                and trial_dist == 0):
            happy = True
        elif gc_stats['gc_cont'] > MAX_GC:
            log("GC content still too high\n")

    log("\nTarget conditions met.\n")
    log(f"Total attempts: {attempts}\n{best_sol_seq}\n")

    _, kl_repeats_val = pattern_prevent(trial_sol, target_nopk, constraint, strand_length)
    log(f"\n{kl_repeats_val} KL repeats found.\n\n")

    # ============================================================
    # PHASE 4: KL Repeat Removal
    # ============================================================

    log("\n__________________________________________________________\n")
    log("Removing KL Repeats\n")

    duplex_cache: dict = {}
    kl_analysis: dict = {}
    best_kl_score = 1000
    kl_score = 0
    happy = False
    pattern_zones = ['-'] * strand_length

    while not happy:
        keep_sequence = False
        attempts += 1
        log(f"({attempts})")
        old_kl_repeats = kl_repeats_val

        for i in range(strand_length):
            if target[i] in ('[', ']') and constraint[i] == 'N' and pattern_zones[i] != '-':
                rand1 = random.uniform(0, 100)
                bp_dir = random.randint(0, 1)
                bp_type = random.uniform(0, 100)
                if kl_repeats_val > 0 and rand1 < 100 * (1 / (4 * kl_repeats_val)):
                    log("!")
                    mi = pair_map[i]
                    if bp_type < KL_GC:
                        trial_sol[i], trial_sol[mi] = ('G', 'C') if bp_dir == 0 else ('C', 'G')
                    else:
                        trial_sol[i], trial_sol[mi] = ('A', 'U') if bp_dir == 0 else ('U', 'A')

        new_trial_seq = ''.join(trial_sol)
        new_fold_str, _ = cached_fold(new_trial_seq)
        trial_dist = rna_distance(new_fold_str, target_str)

        pattern_zones, kl_repeats_val = pattern_prevent(
            trial_sol, target_nopk, constraint, strand_length)
        rpt2 = count_repeats(trial_sol, target, constraint, strand_length,
                             COMPLEMENT_WINDOW, DUPLICATE_WINDOW, PALINDROME_WINDOW)
        pr_cnt = rpt2['pattern_repeats']
        poly_cnt = rpt2['poly_repeats']
        rs_cnt = rpt2['restriction_sites']
        cz_cnt = rpt2['complement_zones']
        dz_cnt = rpt2['duplication_zones']
        pal_cnt = rpt2['palindromes']

        gs_trial = count_gc(trial_sol, target, init_target, pair_map, strand_length)
        if (trial_dist == 0 and kl_repeats_val <= old_kl_repeats and
                gs_trial['gc_cont'] < MAX_GC and pr_cnt == 0 and poly_cnt == 0 and
                rs_cnt == 0 and cz_cnt == 0 and dz_cnt == 0 and pal_cnt == 0):
            keep_sequence = True

        if keep_sequence:
            for i in range(strand_length):
                mut_map[i] = '-' if best_sol[i] == trial_sol[i] else '*'
                best_sol[i] = trial_sol[i]

            best_sol_seq = ''.join(best_sol)
            parent_dist = trial_dist
            parent_fold_str = new_fold_str
            parent_fold = list(parent_fold_str)

            log(f"\nIteration# {attempts}\n")
            gc_stats = log_gc(best_sol)
            log(f"KL Repeats: {kl_repeats_val}\n")
            log(f"{best_sol_seq}\n\n")
        else:
            log(" X \n")
            trial_sol = list(best_sol)

        if trial_dist == 0 and kl_repeats_val == 0:
            happy = True

    log(f"\nInitial design of {name} is complete!\n\n")

    kl_analysis = analyze_kl(
        best_sol, target, target_nopk, pair_map, constraint, strand_length,
        duplex_cache, report_kl=True, full_kl_output=False,
        kl_min_delta_g=KL_MIN_DELTA_G, min_kl=MIN_KL, max_kl=MAX_KL,
        log=log,
    )
    parent_kl_opt_target = list(kl_analysis['kl_opt_target'])

    log("\n\n__________________________________________________________\n")
    log("Optimizing KL Sequences\n")

    happy = False
    best_kl_score = 1000
    kl_score = kl_analysis['kl_score']
    num_kl_sites = kl_analysis['num_kl_sites']

    while not happy:
        attempts += 1
        log(f"\n(Iteration {attempts})")
        log(f" [M{num_kl_sites}] ")

        for i in range(strand_length):
            if target[i] in ('[', ']') and constraint[i] == 'N':
                bp_dir = random.randint(0, 1)
                bp_type = random.uniform(0, 100)
                if random.uniform(0, 1) < (6 * parent_kl_opt_target[i] / (num_kl_sites + 0.01)):
                    log(f"{i}: ")
                    mi = pair_map[i]
                    if bp_type < KL_GC:
                        trial_sol[i], trial_sol[mi] = ('G', 'C') if bp_dir == 0 else ('C', 'G')
                    else:
                        trial_sol[i], trial_sol[mi] = ('A', 'U') if bp_dir == 0 else ('U', 'A')

        new_trial_seq = ''.join(trial_sol)
        pattern_zones, kl_repeats_val = pattern_prevent(
            trial_sol, target_nopk, constraint, strand_length)
        rpt2 = count_repeats(trial_sol, target, constraint, strand_length,
                             COMPLEMENT_WINDOW, DUPLICATE_WINDOW, PALINDROME_WINDOW)
        pr_cnt = rpt2['pattern_repeats']
        poly_cnt = rpt2['poly_repeats']
        rs_cnt = rpt2['restriction_sites']
        cz_cnt = rpt2['complement_zones']
        dz_cnt = rpt2['duplication_zones']
        pal_cnt = rpt2['palindromes']
        repeat_map = rpt2['repeat_map']

        unmasked, masked = _is_masked(repeat_map)
        if unmasked == 0 and masked > 0:
            poly_cnt = 0
            rs_cnt = 0
            pr_cnt = 0
            cz_cnt = 0
            dz_cnt = 0
            pal_cnt = 0

        pattern_prevent_failure = False
        previously_tried = False
        trial_dist = 1
        kl_analyzed = False
        new_fold_str = " "

        if (kl_repeats_val == 0 and poly_cnt == 0 and rs_cnt == 0 and pr_cnt == 0 and
                cz_cnt == 0 and dz_cnt == 0 and pal_cnt == 0):
            log("\nFolding...\n")
            new_fold_str, previously_tried = cached_fold(new_trial_seq)
            trial_dist = rna_distance(new_fold_str, target_str)
        else:
            pattern_prevent_failure = True

        if trial_dist == 0 and kl_repeats_val == 0 and not previously_tried:
            log("\nFolding successful.\nAnalyzing KLs...\n\n")
            kl_analysis = analyze_kl(
                trial_sol, target, target_nopk, pair_map, constraint, strand_length,
                duplex_cache, report_kl=True, full_kl_output=False,
                kl_min_delta_g=KL_MIN_DELTA_G, min_kl=MIN_KL, max_kl=MAX_KL,
                log=log,
            )
            kl_score = kl_analysis['kl_score']
            kl_analyzed = True
        else:
            kl_score = 0

        if (kl_score <= best_kl_score and trial_dist == 0 and kl_repeats_val == 0 and
                not previously_tried and poly_cnt == 0 and rs_cnt == 0 and pr_cnt == 0 and
                cz_cnt == 0 and dz_cnt == 0):
            best_kl_score = kl_score
            parent_kl_opt_target = list(kl_analysis['kl_opt_target'])
            num_kl_sites = kl_analysis['num_kl_sites']

            for i in range(strand_length):
                mut_map[i] = '-' if best_sol[i] == trial_sol[i] else trial_sol[i]
                best_sol[i] = trial_sol[i]

            best_sol_seq = ''.join(best_sol)
            parent_dist = trial_dist
            parent_fold_str = new_fold_str
            parent_fold = list(parent_fold_str)

            log(f"\nIteration# {attempts}\n")
            gc_stats = log_gc(best_sol)
            log(f"KL Repeats: {kl_repeats_val}\n")
            log(f"KL Score: {kl_score}\n")
            log(f"\nKL Target Mask:\n{''.join(str(x) for x in parent_kl_opt_target)}\n\n")
            log(f"{best_sol_seq}\n\n")
            log(f"{''.join(mut_map)}\n\n")
        else:
            log(f"\nRound {attempts} Failed:   ")
            if kl_repeats_val > 0:
                log(f"   {kl_repeats_val} KL repeats were found\n")
            if previously_tried:
                log("   The Sequence was previously tested\n")
            if kl_score > best_kl_score and kl_analyzed:
                log(f"   The KL score {kl_score} is greater than the best score {best_kl_score}\n")
            if pattern_prevent_failure:
                log("   Pattern Prevent Failed.\n")

            if kl_analysis:
                kl_analysis['kl_opt_target'] = list(parent_kl_opt_target)
            trial_sol = list(best_sol)

        if trial_dist == 0 and kl_score == 0:
            happy = True

    print_victory(log, start_time, failed_initial_design)
    log(f"\nOptimization of design {name} is complete!\n\n")
    log(f"\n\n\nIteration# {attempts}\n")
    gc_stats = log_gc(best_sol)
    log(f"KL Repeats: {kl_repeats_val}\nKL Score: {kl_score}\n{best_sol_seq}\n\n")

    analyze_kl(
        trial_sol, target, target_nopk, pair_map, constraint, strand_length,
        duplex_cache, report_kl=False, full_kl_output=True,
        kl_min_delta_g=KL_MIN_DELTA_G, min_kl=MIN_KL, max_kl=MAX_KL,
        log=log,
    )

    # ============================================================
    # PHASE 5: Ensemble Diversity Optimisation
    # ============================================================

    log("\nBeginning Ensemble Defect Optimization.\n\n")
    start_time = time.time()  # reset timer for ED phase

    lowest_ed = 100_000_000.0
    happy     = False
    test_first = True

    rpt = count_repeats(trial_sol, target, constraint, strand_length,
                        COMPLEMENT_WINDOW, DUPLICATE_WINDOW, PALINDROME_WINDOW)
    last_pr  = rpt['pattern_repeats']
    last_poly = rpt['poly_repeats']
    last_rs  = rpt['restriction_sites']
    last_cz  = rpt['complement_zones']
    last_dz  = rpt['duplication_zones']
    last_pal = rpt['palindromes']

    while not happy:
        rpt = count_repeats(trial_sol, target, constraint, strand_length,
                            COMPLEMENT_WINDOW, DUPLICATE_WINDOW, PALINDROME_WINDOW)
        repeat_map = rpt['repeat_map']
        pr_cnt   = rpt['pattern_repeats']
        poly_cnt = rpt['poly_repeats']
        rs_cnt   = rpt['restriction_sites']
        cz_cnt   = rpt['complement_zones']
        dz_cnt   = rpt['duplication_zones']
        pal_cnt  = rpt['palindromes']

        GC_check = 100
        keep_sequence = False

        if not test_first:
            attempts += 1

            while True:
                log(f"(Iteration {attempts})")
                num_mutation_spots = 0
                GC_check = 100
                mutation_sites = [0] * strand_length

                for i in range(strand_length):
                    trial_sol[i] = best_sol[i]

                gs = count_gc(best_sol, target, init_target, pair_map, strand_length)
                gc_pct = gs['gc_pairs'] / max(1, gs['pairs']) * 100

                for i in range(strand_length):
                    rand1 = random.uniform(0, 100)
                    if target[i] in ('(', '.') and rand1 < drift_rate * 4:
                        mutation_sites[i] += 1
                        num_mutation_spots += 1

                if rad_level < 2:
                    rad_level = 2
                if rad_level > 20:
                    rad_level = 20

                mut_scale = rad_level / (num_mutation_spots + 0.001)
                log(f" M[{num_mutation_spots}] (*{rad_level}) ")

                for i in range(strand_length - 1):
                    if random.uniform(0, 1) >= mut_scale * mutation_sites[i]:
                        continue

                    ci = constraint[i]
                    mi = pair_map[i]
                    cmi = constraint[mi] if mi != i else 'N'
                    ti = target[i]
                    coin = random.randint(0, 1)
                    a, b = trial_sol[i], trial_sol[mi]
                    rm = repeat_map[i]

                    if ti == '.' and ci == 'N':
                        log('.')
                        trial_sol[i] = rand_base()

                    if (rm in ('G', 'A') and ci == 'R' and random.uniform(0, 1) < 0.5 and
                            ti in ('(', ')', '[', ']') and cmi in ('N', 'Y')):
                        log(':R ')
                        if a == 'G' and (gc_pct > MAX_GC or random.uniform(0, 1) < 0.5):
                            trial_sol[i] = 'A'
                            trial_sol[mi] = 'U'
                        else:
                            trial_sol[i] = 'G'
                            trial_sol[mi] = 'C'

                    if (rm in ('X', 'U', 'C') and ci == 'Y' and random.uniform(0, 1) < 0.5 and
                            ti in ('(', ')', '[', ']') and cmi in ('N', 'R')):
                        log(':Y ')
                        if a == 'U' and (gc_pct > MAX_GC or random.uniform(0, 1) < 0.5):
                            trial_sol[i] = 'C'
                            trial_sol[mi] = 'G'
                        else:
                            trial_sol[i] = 'U'
                            trial_sol[mi] = 'A'

                    if ti == '(' and ci == 'N' and cmi == 'N':
                        log(':')
                        if {a, b} == {'C', 'G'}:
                            if gc_pct > MAX_GC or random.uniform(0, 1) < 0.5:
                                log('+A ')
                                trial_sol[i], trial_sol[mi] = ('A', 'U') if coin == 0 else ('U', 'A')
                        elif {a, b} == {'A', 'U'}:
                            if gc_pct < MAX_GC or random.uniform(0, 1) < 0.5:
                                log('+G ')
                                trial_sol[i], trial_sol[mi] = ('G', 'C') if coin == 0 else ('C', 'G')
                        elif {a, b} == {'G', 'U'}:
                            if random.uniform(0, 1) < 0.9:
                                log(':fW ')
                                trial_sol[i], trial_sol[mi] = ('U', 'G') if a == 'G' else ('G', 'U')
                            else:
                                log(':xW ')
                                if a == 'G':
                                    trial_sol[i] = 'A'
                                else:
                                    trial_sol[mi] = 'A'
                    elif ti in ('(', ')') and ci == 'K' and cmi == 'K':
                        log(':fW ')
                        trial_sol[i], trial_sol[mi] = ('U', 'G') if a == 'G' else ('G', 'U')
                    elif ti in ('(', ')') and ci == 'S' and cmi == 'S':
                        log(':fS ')
                        trial_sol[i], trial_sol[mi] = ('C', 'G') if a == 'G' else ('G', 'C')

                gs2 = count_gc(trial_sol, target, init_target, pair_map, strand_length)
                GC_check = 100 if gs2['gc_cont'] > MAX_GC_OPT else 0

                rpt2 = count_repeats(trial_sol, target, constraint, strand_length,
                                     COMPLEMENT_WINDOW, DUPLICATE_WINDOW, PALINDROME_WINDOW)
                pr_cnt = rpt2['pattern_repeats']
                poly_cnt = rpt2['poly_repeats']
                rs_cnt = rpt2['restriction_sites']
                cz_cnt = rpt2['complement_zones']
                dz_cnt = rpt2['duplication_zones']
                pal_cnt = rpt2['palindromes']
                repeat_map = rpt2['repeat_map']

                log(f"(Repeats:{pr_cnt}/Poly:{poly_cnt}/Restrict:{rs_cnt}/Needs Wobble:{cz_cnt}/GC Content:{gs2['gc_cont']} Palindromes:{pal_cnt})\n")

                unmasked, masked = _is_masked(repeat_map)
                if unmasked == 0 and masked > 0:
                    log("\n\n Remaining patterns are masked by sequence constraints.  Passing to next phase. \n\n")
                    pr_cnt = poly_cnt = rs_cnt = cz_cnt = dz_cnt = pal_cnt = 0
                    failed_pattern_test = True

                cur = _pattern_score(pr_cnt, poly_cnt, rs_cnt, cz_cnt, dz_cnt, pal_cnt)
                lst = _pattern_score(last_pr, last_poly, last_rs, last_cz, last_dz, last_pal)
                if cur > lst:
                    rad_level -= 1
                if GC_check > 0:
                    rad_level += 6

                if pr_cnt + poly_cnt + rs_cnt + cz_cnt + dz_cnt + pal_cnt + GC_check == 0:
                    break

        test_first = False

        new_trial_seq = ''.join(trial_sol)
        log("SUCCESS folding...")
        new_fold_str, _ = cached_fold(new_trial_seq)
        trial_dist = rna_distance(new_fold_str, target_str)

        gs_trial = count_gc(trial_sol, target, init_target, pair_map, strand_length)
        GC_check = 100 if gs_trial['gc_cont'] > MAX_GC_OPT else 0

        if (trial_dist == 0 and GC_check == 0 and
                (pr_cnt <= last_pr or poly_cnt <= last_poly or rs_cnt <= last_rs or
                 cz_cnt > last_cz or dz_cnt > last_dz or pal_cnt > last_pal)):
            keep_sequence = True

        keep_ed = False
        if keep_sequence:
            log("\nKeeping sequence\n")
            lowest_ed, keep_ed = export_result(
                trial_sol, name, output_file_path, target_pk_str,
                gs_trial, kl_analysis, failed_initial_design, failed_pattern_test,
                lowest_ed, log,
            )
            if lowest_ed <= TARGET_ED:
                log(f" Ensemble Diversity Target Reached. {lowest_ed} <= {TARGET_ED}\n")
                happy = True
        else:
            log(" Mutant does not meet design criteria...\n\n")

        if keep_sequence and keep_ed:
            log("Seeding next generation with lowest ensemble diversity design.\n\n")
            last_pr  = pr_cnt; last_poly = poly_cnt; last_rs = rs_cnt
            last_cz  = cz_cnt; last_dz   = dz_cnt;  last_pal = pal_cnt
            rad_level += 3

            for i in range(strand_length):
                mut_map[i]  = trial_sol[i] if best_sol[i] != trial_sol[i] else '-'
                best_sol[i] = trial_sol[i]

            best_sol_seq    = ''.join(best_sol)
            parent_dist     = trial_dist
            parent_fold_str = new_fold_str
            parent_fold     = list(parent_fold_str)
        else:
            rad_level -= 1
            if rad_level < 3:
                rad_level = 20  # escape local minimum by jumping to high radiation
            if trial_dist > 0:
                log(f" Bad Fold: -D{trial_dist} ")
            log("Reverting to previous design.\n\n")
            trial_sol = list(best_sol)

    spool_file.close()
    print(f"\nDone. Full log written to {spool_path}")


if __name__ == '__main__':
    main()
