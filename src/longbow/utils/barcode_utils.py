import gzip
import sys
import enum


import symspellpy
from symspellpy import SymSpell, Verbosity

from tqdm import tqdm

from ordered_set import OrderedSet


class SymSpellMatchResultType(enum.Enum):
    SUCCESS = enum.auto()
    SUCCESS_BARCODE_ALREADY_CORRECT = enum.auto()
    AMBIGUOUS = enum.auto()
    NO_MATCH_IN_LEV_DIST = enum.auto()


def _load_barcode_allowlist_helper(open_file_handle, barcode_set, disable_pbar=None):
    if disable_pbar is None:
        disable_pbar = not sys.stdin.isatty()
    for line in tqdm(open_file_handle,
                     desc="Loading barcode allow list", unit=" barcode", colour="green", file=sys.stderr, leave=True,
                     disable=disable_pbar):
        bc = line.decode("utf-8").strip()
        barcode_set.add(bc)


def load_barcode_allowlist(filename, disable_pbar=None):
    barcodes = OrderedSet()

    if filename is not None:
        if filename.endswith(".gz"):
            with gzip.open(filename, 'rb') as f:
                _load_barcode_allowlist_helper(f, barcodes, disable_pbar)
        else:
            with open(filename, 'rb') as f:
                _load_barcode_allowlist_helper(f, barcodes, disable_pbar)

    return barcodes


def generate_symspell_index(wl_fname, max_dist_allowed, prefix_len=16, freq_list_sep="\t"):
    # Code generously donated from Victoria Popic
    sym_spell = SymSpell(max_dictionary_edit_distance=max_dist_allowed, prefix_length=prefix_len)
    sym_spell._distance_algorithm = symspellpy.editdistance.DistanceAlgorithm.LEVENSHTEIN_FAST
    sym_spell.load_dictionary(wl_fname, 0, 1, separator=freq_list_sep)  # positions of the term and count, respectively
    return sym_spell


def find_match_symspell(barcode, barcode_allow_list, sym_spell_index, dist_thr):
    # Adapted from code generously donated from Victoria Popic
    # using faster set lookups to find exact matches
    # for the most common case
    if barcode in barcode_allow_list:  # there is an exact match
        return barcode, 0, SymSpellMatchResultType.SUCCESS_BARCODE_ALREADY_CORRECT

    matches = sym_spell_index.lookup(barcode, Verbosity.CLOSEST, max_edit_distance=dist_thr)
    if not matches:
        return None, None, SymSpellMatchResultType.NO_MATCH_IN_LEV_DIST
    elif len(matches) > 1:
        return None, None, SymSpellMatchResultType.AMBIGUOUS

    return matches[0].term, matches[0].distance, SymSpellMatchResultType.SUCCESS
