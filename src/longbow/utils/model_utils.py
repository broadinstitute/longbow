import sys
import re
import logging

import click_log

from pomegranate import State
from pomegranate import *

import longbow.utils.bam_utils
import longbow.utils.constants
from .constants import RANDOM_SEGMENT_NAME, FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME, HPR_SEGMENT_TYPE_NAME, \
    RANDOM_SILENT_STATE_A, RANDOM_SILENT_STATE_B, RANDOM_BASE_STATE, START_STATE_INDICATOR, END_STATE_INDICATOR, \
    BAKE_MERGE_STRATEGY

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger(__name__)
click_log.basic_config(logger)

starts_with_number_re = re.compile(r"^\d")


class ModelBuilder:
    """Utilities for constructing a full Longbow model."""

    # Define constants for all our default probabilities here:
    RANDOM_BASE_PROB = 0.25

    PER_BASE_MATCH_PROB = 0.94
    PER_BASE_MISMATCH_PROB = 0.02

    MATCH_MATCH_PROB = 0.90
    MATCH_INDEL_PROB = 0.05
    MATCH_TRAIL_INSERT_PROB = 0.90

    INDEL_CONTINUATION_PROB = 0.7
    INDEL_SWITCH_PROB = 0.3

    START_RANDOM_PROB = 1.0

    START_AND_END_RANDOM_PROB = 0.5
    RAND_RAND_PROB = 0.5
    RAND_INS_CONTINUATION_PROB = 0.8
    RAND_INS_TO_DEL_PROB = 0.1
    RAND_INS_END_PROB = 0.1

    NAMED_RAND_CONTINUE_PROB = 0.9
    NAMED_RAND_EXIT_PROB = 1-NAMED_RAND_CONTINUE_PROB

    HPR_MATCH_PROB = 0.9
    HPR_MISMATCH_PROB = (1-HPR_MATCH_PROB) / 3

    HPR_BOOKEND_MATCH_PROB = 0.99
    HPR_BOOKEND_MISMATCH_PROB = (1-HPR_BOOKEND_MATCH_PROB)/3

    HPR_SUDDEN_END_PROB = 0.01
    HPR_MODEL_RECURRENCE_PROB = 0.1

    SUDDEN_END_PROB = 0.01
    MATCH_END_PROB = 0.1


    @staticmethod
    def make_global_alignment_model(target, name=None):
        logger.debug("Making Model: GLOBAL_ALIGNMENT (%s)", name)
        model = HiddenMarkovModel(name=name)
        s = {}

        # add states
        i0 = State(
            DiscreteDistribution({
                "A": ModelBuilder.RANDOM_BASE_PROB,
                "C": ModelBuilder.RANDOM_BASE_PROB,
                "G": ModelBuilder.RANDOM_BASE_PROB,
                "T": ModelBuilder.RANDOM_BASE_PROB
            }),
            name=f"{name}:I0",
        )

        model.add_state(i0)

        s[i0.name] = i0

        for c in range(len(target)):
            dc = State(None, name=f"{name}:D{c + 1}")

            mc = State(
                DiscreteDistribution({
                    "A": ModelBuilder.PER_BASE_MATCH_PROB if target[c] == "A" else ModelBuilder.PER_BASE_MISMATCH_PROB,
                    "C": ModelBuilder.PER_BASE_MATCH_PROB if target[c] == "C" else ModelBuilder.PER_BASE_MISMATCH_PROB,
                    "G": ModelBuilder.PER_BASE_MATCH_PROB if target[c] == "G" else ModelBuilder.PER_BASE_MISMATCH_PROB,
                    "T": ModelBuilder.PER_BASE_MATCH_PROB if target[c] == "T" else ModelBuilder.PER_BASE_MISMATCH_PROB,
                }),
                name=f"{name}:M{c + 1}",
            )

            ic = State(
                DiscreteDistribution({
                    "A": ModelBuilder.RANDOM_BASE_PROB,
                    "C": ModelBuilder.RANDOM_BASE_PROB,
                    "G": ModelBuilder.RANDOM_BASE_PROB,
                    "T": ModelBuilder.RANDOM_BASE_PROB
                }),
                name=f"{name}:I{c + 1}",
            )

            model.add_states([mc, ic, dc])

            s[dc.name] = dc
            s[mc.name] = mc
            s[ic.name] = ic

        # add transitions
        model.add_transition(model.start, s[f"{name}:I0"], ModelBuilder.MATCH_INDEL_PROB)
        model.add_transition(model.start, s[f"{name}:D1"], ModelBuilder.MATCH_INDEL_PROB)
        model.add_transition(model.start, s[f"{name}:M1"], ModelBuilder.MATCH_MATCH_PROB)

        model.add_transition(s[f"{name}:I0"], s[f"{name}:I0"], ModelBuilder.INDEL_CONTINUATION_PROB)
        model.add_transition(s[f"{name}:I0"], s[f"{name}:D1"], ModelBuilder.INDEL_SWITCH_PROB/2)
        model.add_transition(s[f"{name}:I0"], s[f"{name}:M1"], ModelBuilder.INDEL_SWITCH_PROB/2)

        for c in range(1, len(target)):
            model.add_transition(s[f"{name}:D{c}"], s[f"{name}:D{c + 1}"], ModelBuilder.INDEL_SWITCH_PROB/2)
            model.add_transition(s[f"{name}:D{c}"], s[f"{name}:I{c}"], ModelBuilder.INDEL_CONTINUATION_PROB)
            model.add_transition(s[f"{name}:D{c}"], s[f"{name}:M{c + 1}"], ModelBuilder.INDEL_SWITCH_PROB/2)

            model.add_transition(s[f"{name}:I{c}"], s[f"{name}:D{c + 1}"], ModelBuilder.INDEL_SWITCH_PROB/2)
            model.add_transition(s[f"{name}:I{c}"], s[f"{name}:I{c}"], ModelBuilder.INDEL_SWITCH_PROB/2)
            model.add_transition(s[f"{name}:I{c}"], s[f"{name}:M{c + 1}"], ModelBuilder.INDEL_CONTINUATION_PROB)

            model.add_transition(s[f"{name}:M{c}"], s[f"{name}:D{c + 1}"], ModelBuilder.MATCH_INDEL_PROB)
            model.add_transition(s[f"{name}:M{c}"], s[f"{name}:I{c}"], ModelBuilder.MATCH_INDEL_PROB)
            model.add_transition(s[f"{name}:M{c}"], s[f"{name}:M{c + 1}"], ModelBuilder.MATCH_MATCH_PROB)

        model.add_transition(
            s[f"{name}:D{len(target)}"], s[f"{name}:I{len(target)}"], ModelBuilder.INDEL_CONTINUATION_PROB
        )
        model.add_transition(s[f"{name}:D{len(target)}"], model.end, ModelBuilder.INDEL_SWITCH_PROB)

        model.add_transition(s[f"{name}:I{len(target)}"], s[f"{name}:I{len(target)}"], ModelBuilder.INDEL_SWITCH_PROB/2)
        model.add_transition(
            s[f"{name}:I{len(target)}"],
            model.end,
            ModelBuilder.INDEL_CONTINUATION_PROB + (ModelBuilder.INDEL_SWITCH_PROB/2)
        )

        model.add_transition(
            s[f"{name}:M{len(target)}"], s[f"{name}:I{len(target)}"], ModelBuilder.MATCH_TRAIL_INSERT_PROB
        )
        model.add_transition(s[f"{name}:M{len(target)}"], model.end, ModelBuilder.MATCH_END_PROB)

        model.bake(merge=BAKE_MERGE_STRATEGY)
        return model

    @staticmethod
    def make_homopolymer_repeat_model(name, nucleotide, expected_length):
        logger.debug("Making Model: HOMOPOLYMER_REPEAT (%s:%s x %d)", name, nucleotide, expected_length)

        model = HiddenMarkovModel(name=name)
        s = {}

        # Define our distributions here so we can use them later:
        random_distribution = DiscreteDistribution({
            "A": ModelBuilder.RANDOM_BASE_PROB,
            "C": ModelBuilder.RANDOM_BASE_PROB,
            "G": ModelBuilder.RANDOM_BASE_PROB,
            "T": ModelBuilder.RANDOM_BASE_PROB
        })

        standard_distribution = DiscreteDistribution({
            "A": ModelBuilder.PER_BASE_MATCH_PROB if nucleotide == "A" else ModelBuilder.PER_BASE_MISMATCH_PROB,
            "C": ModelBuilder.PER_BASE_MATCH_PROB if nucleotide == "C" else ModelBuilder.PER_BASE_MISMATCH_PROB,
            "G": ModelBuilder.PER_BASE_MATCH_PROB if nucleotide == "G" else ModelBuilder.PER_BASE_MISMATCH_PROB,
            "T": ModelBuilder.PER_BASE_MATCH_PROB if nucleotide == "T" else ModelBuilder.PER_BASE_MISMATCH_PROB,
        })

        # Add an extra penalty to the last state for being a base that isn't the one in this HPR:
        bookend_state_distribution = DiscreteDistribution({
            "A": ModelBuilder.HPR_BOOKEND_MATCH_PROB if nucleotide == "A" else ModelBuilder.HPR_BOOKEND_MISMATCH_PROB,
            "C": ModelBuilder.HPR_BOOKEND_MATCH_PROB if nucleotide == "C" else ModelBuilder.HPR_BOOKEND_MISMATCH_PROB,
            "G": ModelBuilder.HPR_BOOKEND_MATCH_PROB if nucleotide == "G" else ModelBuilder.HPR_BOOKEND_MISMATCH_PROB,
            "T": ModelBuilder.HPR_BOOKEND_MATCH_PROB if nucleotide == "T" else ModelBuilder.HPR_BOOKEND_MISMATCH_PROB,
        })

        # add states
        i0 = State(random_distribution, name=f"{name}:I0")

        model.add_state(i0)
        s[i0.name] = i0

        for c in range(expected_length):
            if c == 0 or c == expected_length-1:
                mc = State(bookend_state_distribution, name=f"{name}:M{c + 1}")
            else:
                mc = State(standard_distribution, name=f"{name}:M{c + 1}")

            ic = State(random_distribution, name=f"{name}:I{c + 1}")

            model.add_states([mc, ic])

            s[mc.name] = mc
            s[ic.name] = ic

        # Add transitions for starting states:
        model.add_transition(model.start, s[f"{name}:I0"], ModelBuilder.MATCH_INDEL_PROB)
        model.add_transition(model.start, s[f"{name}:M1"], ModelBuilder.MATCH_MATCH_PROB)
        model.add_transition(s[f"{name}:M1"], model.end, ModelBuilder.HPR_SUDDEN_END_PROB)
        model.add_transition(s[f"{name}:M1"], model.start, ModelBuilder.HPR_MODEL_RECURRENCE_PROB)

        model.add_transition(s[f"{name}:I0"], s[f"{name}:I0"], ModelBuilder.INDEL_CONTINUATION_PROB)
        model.add_transition(s[f"{name}:I0"], s[f"{name}:M1"], ModelBuilder.INDEL_SWITCH_PROB / 2)

        # Add transitions for middle states:
        for c in range(1, expected_length):
            model.add_transition(s[f"{name}:I{c}"], s[f"{name}:I{c}"], ModelBuilder.INDEL_SWITCH_PROB / 2)
            model.add_transition(s[f"{name}:I{c}"], s[f"{name}:M{c + 1}"], ModelBuilder.INDEL_CONTINUATION_PROB)
            model.add_transition(s[f"{name}:I{c}"], model.start, ModelBuilder.HPR_MODEL_RECURRENCE_PROB)
            model.add_transition(s[f"{name}:I{c}"], model.end, ModelBuilder.HPR_SUDDEN_END_PROB)

            model.add_transition(s[f"{name}:M{c}"], s[f"{name}:I{c}"], ModelBuilder.MATCH_INDEL_PROB)
            model.add_transition(s[f"{name}:M{c}"], s[f"{name}:M{c + 1}"], ModelBuilder.MATCH_MATCH_PROB)

        # Add transitions for the last states:
        model.add_transition(s[f"{name}:I{expected_length}"], s[f"{name}:I{expected_length}"], ModelBuilder.INDEL_SWITCH_PROB / 2)
        model.add_transition(
            s[f"{name}:I{expected_length}"],
            model.end,
            ModelBuilder.INDEL_CONTINUATION_PROB + (ModelBuilder.INDEL_SWITCH_PROB / 2)
        )

        model.add_transition(
            s[f"{name}:M{expected_length}"], s[f"{name}:I{expected_length}"], ModelBuilder.MATCH_TRAIL_INSERT_PROB
        )
        model.add_transition(s[f"{name}:M{expected_length}"], model.end, ModelBuilder.MATCH_END_PROB)
        model.add_transition(s[f"{name}:M{expected_length}"], model.start, ModelBuilder.MATCH_END_PROB)

        # Finalize and return the model:
        model.bake(merge=BAKE_MERGE_STRATEGY)
        return model

    @staticmethod
    def make_named_random_model(name):
        logger.debug("Making Model: NAMED RANDOM (%s)", name)

        model = HiddenMarkovModel(name=name)

        # This is a VERY simple repeat model.
        # The idea here is that is a user specifies a named random segment, then we know
        # there should be at least one random base here.
        ri = State(
            DiscreteDistribution({
                "A": ModelBuilder.RANDOM_BASE_PROB,
                "C": ModelBuilder.RANDOM_BASE_PROB,
                "G": ModelBuilder.RANDOM_BASE_PROB,
                "T": ModelBuilder.RANDOM_BASE_PROB
            }),
            name=f"{name}:{RANDOM_BASE_STATE}",
        )

        model.add_state(ri)

        model.add_transition(model.start, ri, 1)
        model.add_transition(ri, ri, ModelBuilder.NAMED_RAND_CONTINUE_PROB)
        model.add_transition(ri, model.end, ModelBuilder.NAMED_RAND_EXIT_PROB)

        model.bake(merge=BAKE_MERGE_STRATEGY)
        return model

    @staticmethod
    def make_random_repeat_model():
        logger.debug(f"Making Model: RANDOM REPEAT")
        model = HiddenMarkovModel(name=RANDOM_SEGMENT_NAME)

        # add states
        ri = State(
            DiscreteDistribution({
                "A": ModelBuilder.RANDOM_BASE_PROB,
                "C": ModelBuilder.RANDOM_BASE_PROB,
                "G": ModelBuilder.RANDOM_BASE_PROB,
                "T": ModelBuilder.RANDOM_BASE_PROB
            }),
            name=f"{RANDOM_SEGMENT_NAME}:{RANDOM_BASE_STATE}",
        )
        rda = State(None, name=f"{RANDOM_SEGMENT_NAME}:{RANDOM_SILENT_STATE_A}")
        rdb = State(None, name=f"{RANDOM_SEGMENT_NAME}:{RANDOM_SILENT_STATE_B}")

        model.add_states([ri, rda, rdb])

        # add transitions
        model.add_transition(model.start, rda, ModelBuilder.START_AND_END_RANDOM_PROB)
        model.add_transition(model.start, ri, ModelBuilder.START_AND_END_RANDOM_PROB)

        model.add_transition(ri, ri, ModelBuilder.RAND_INS_CONTINUATION_PROB)
        model.add_transition(ri, rda, ModelBuilder.RAND_INS_TO_DEL_PROB)
        model.add_transition(ri, model.end, ModelBuilder.RAND_INS_END_PROB)

        model.add_transition(rdb, ri, ModelBuilder.RAND_RAND_PROB)
        model.add_transition(rdb, model.end, ModelBuilder.START_AND_END_RANDOM_PROB)

        model.bake(merge=BAKE_MERGE_STRATEGY)

        return model

    @staticmethod
    def make_fixed_length_random_segment(name, length):
        logger.debug("Making Model: FIXED_LENGTH_RANDOM (%s:%d)", name, length)
        model = HiddenMarkovModel(name=name)

        state_dict = {}

        random_dist = DiscreteDistribution({
            "A": ModelBuilder.RANDOM_BASE_PROB,
            "C": ModelBuilder.RANDOM_BASE_PROB,
            "G": ModelBuilder.RANDOM_BASE_PROB,
            "T": ModelBuilder.RANDOM_BASE_PROB
        })

        # Create match states for each base in the fixed length random segment:
        for i in range(length):
            match_state = State(random_dist, name=f"{name}:M{i + 1}")

            model.add_state(match_state)
            state_dict[match_state.name] = match_state

        # Add transitions:
        # NOTE: All probabilities in the transitions are linear because we specifically want to capture a chain of
        # fixed-length here, so each state as exactly one transition.
        model.add_transition(model.start, state_dict[f"{name}:M1"], 1)
        for i in range(1, length):
            model.add_transition(state_dict[f"{name}:M{i}"], state_dict[f"{name}:M{i + 1}"], 1)
        model.add_transition(state_dict[f"{name}:M{length}"], model.end, 1)

        # Finalize and return the model:
        model.bake(merge=BAKE_MERGE_STRATEGY)
        return model

    @staticmethod
    def find_state(base_hmm, state_name):
        for state in base_hmm.states:
            if state.name == state_name:
                return state

        return None

    @staticmethod
    def connect_terminals(base_hmm, adapter_name_i, adapter_name_j, transition_probability=1.0):
        prev_state = ModelBuilder.find_state(base_hmm, f'{adapter_name_i}-end')
        next_state = ModelBuilder.find_state(base_hmm, f'{adapter_name_j}-start')

        if prev_state is not None and next_state is not None:
            base_hmm.add_transition(prev_state, next_state, transition_probability)

        # base_hmm.bake(merge=BAKE_MERGE_STRATEGY)
        return base_hmm

    pre_configured_models = {
        'array': {
            "mas_10": {
                "description": "10-element MAS-ISO-seq array",
                "version": "3.0.0",
                "structure": [ "Q", "C", "M", "I", "O", "J", "B", "D", "K", "H", "R" ],
                "adapters": {
                    "Q": "AAGCACCATAATGTGT",
                    "C": "ACTCTGTCAGGTCCGA",
                    "M": "ACCTAGATCAGAGCCT",
                    "I": "AGTGCGTTGCGAATTG",
                    "O": "AAGTCACCGGCACCTT",
                    "J": "AATTGCGTAGTTGGCC",
                    "B": "ACTTGTAAGCTGTCTA",
                    "D": "ACCTCCTCCTCCAGAA",
                    "K": "ACACTTGGTCGCAATC",
                    "H": "ATGTTGAATCCTAGCG",
                    "R": "AACCGGACACACTTAG",
                },
                "deprecated": False,
            },

            "mas_15": {
                "description": "15-element MAS-ISO-seq array",
                "version": "3.0.0",
                "structure": [ "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P" ],
                "adapters": {
                    "A": "AGCTTACTTGTGAAGA",
                    "B": "ACTTGTAAGCTGTCTA",
                    "C": "ACTCTGTCAGGTCCGA",
                    "D": "ACCTCCTCCTCCAGAA",
                    "E": "AACCGGACACACTTAG",
                    "F": "AGAGTCCAATTCGCAG",
                    "G": "AATCAAGGCTTAACGG",
                    "H": "ATGTTGAATCCTAGCG",
                    "I": "AGTGCGTTGCGAATTG",
                    "J": "AATTGCGTAGTTGGCC",
                    "K": "ACACTTGGTCGCAATC",
                    "L": "AGTAAGCCTTCGTGTC",
                    "M": "ACCTAGATCAGAGCCT",
                    "N": "AGGTATGCCGGTTAAG",
                    "O": "AAGTCACCGGCACCTT",
                    "P": "ATGAAGTGGCTCGAGA",
                },
                "deprecated": False,
            },

            "mas_16": {
                "description": "16-element MAS-ISO-seq array",
                "version": "3.0.0",
                "structure": [ "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q" ],
                "adapters": {
                    "A": "AGCTTACTTGTGAAGA",
                    "B": "ACTTGTAAGCTGTCTA",
                    "C": "ACTCTGTCAGGTCCGA",
                    "D": "ACCTCCTCCTCCAGAA",
                    "E": "AACCGGACACACTTAG",
                    "F": "AGAGTCCAATTCGCAG",
                    "G": "AATCAAGGCTTAACGG",
                    "H": "ATGTTGAATCCTAGCG",
                    "I": "AGTGCGTTGCGAATTG",
                    "J": "AATTGCGTAGTTGGCC",
                    "K": "ACACTTGGTCGCAATC",
                    "L": "AGTAAGCCTTCGTGTC",
                    "M": "ACCTAGATCAGAGCCT",
                    "N": "AGGTATGCCGGTTAAG",
                    "O": "AAGTCACCGGCACCTT",
                    "P": "ATGAAGTGGCTCGAGA",
                    "Q": "AGTAGCTGTGTGCA",
                },
                "deprecated": False,
            },

            "mas_debug": {
                "description": "All available MAS-ISO-seq adapters in both orientations",
                "version": "3.0.0",
                "structure": [ "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "V",
                               "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "v", ],
                "adapters": {
                    "A": "AGCTTACTTGTGAAGA",
                    "B": "ACTTGTAAGCTGTCTA",
                    "C": "ACTCTGTCAGGTCCGA",
                    "D": "ACCTCCTCCTCCAGAA",
                    "E": "AACCGGACACACTTAG",
                    "F": "AGAGTCCAATTCGCAG",
                    "G": "AATCAAGGCTTAACGG",
                    "H": "ATGTTGAATCCTAGCG",
                    "I": "AGTGCGTTGCGAATTG",
                    "J": "AATTGCGTAGTTGGCC",
                    "K": "ACACTTGGTCGCAATC",
                    "L": "AGTAAGCCTTCGTGTC",
                    "M": "ACCTAGATCAGAGCCT",
                    "N": "AGGTATGCCGGTTAAG",
                    "O": "AAGTCACCGGCACCTT",
                    "P": "ATGAAGTGGCTCGAGA",
                    "Q": "AGTAGCTGTGTGCA",
                    "R": "AACCGGACACACTTAG",
                    "V": "AAGCACCATAATGTGT",

                    "a": "TCTTCACAAGTAAGCT",
                    "b": "TAGACAGCTTACAAGT",
                    "c": "TCGGACCTGACAGAGT",
                    "d": "TTCTGGAGGAGGAGGT",
                    "e": "CTAAGTGTGTCCGGTT",
                    "f": "CTGCGAATTGGACTCT",
                    "g": "CCGTTAAGCCTTGATT",
                    "h": "CGCTAGGATTCAACAT",
                    "i": "CAATTCGCAACGCACT",
                    "j": "GGCCAACTACGCAATT",
                    "k": "GATTGCGACCAAGTGT",
                    "l": "GACACGAAGGCTTACT",
                    "m": "AGGCTCTGATCTAGGT",
                    "n": "CTTAACCGGCATACCT",
                    "o": "AAGGTGCCGGTGACTT",
                    "p": "TCTCGAGCCACTTCAT",
                    "q": "TGCACACAGCTACT",
                    "r": "CTAAGTGTGTCCGGTT",
                    "v": "ACACATTATGGTGCTT",
                },
                "deprecated": False,
            },

            "isoseq": {
                "description": "PacBio IsoSeq model",
                "version": "3.0.0",
                "array_element_structure": [ "V", "M" ],
                "adapters": {
                    "V": "TCTACACGACGCTCTTCCGATCT",
                    "M": "GTACTCTGCGTTGATACCACTGCTT",
                },
                "deprecated": False,
            },
        },

        'cdna': {
            "sc_10x3p": {
                "description": "single-cell 10x 3' kit",
                "version": "3.0.0",
                "structure": [ "5p_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "3p_Adapter" ],
                "adapters": {
                    "5p_Adapter": "TCTACACGACGCTCTTCCGATCT",
                    "CBC": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 16},
                    "UMI": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 12},
                    "Poly_T": {HPR_SEGMENT_TYPE_NAME: ("T", 30)},
                    "cDNA": RANDOM_SEGMENT_NAME,
                    "3p_Adapter": "CCCATGTACTCTGCGTTGATACCACTGCTT",
                },
                "named_random_segments": ["CBC", "UMI", "cDNA"],
                "coding_region": "cDNA",
                "annotation_segments": {
                    "UMI": [(longbow.utils.constants.READ_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG),
                            (longbow.utils.constants.READ_RAW_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG)],
                    "CBC": [(longbow.utils.constants.READ_BARCODE_TAG, longbow.utils.constants.READ_BARCODE_POS_TAG),
                            (longbow.utils.constants.READ_RAW_BARCODE_TAG, longbow.utils.constants.READ_BARCODE_POS_TAG)],
                },
                "deprecated": False,
            },

            "sc_10x5p": {
                "description": "single-cell 10x 5' kit",
                "version": "3.0.0",
                "structure": [ "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter" ],
                "adapters": {
                    "5p_Adapter": "TCTACACGACGCTCTTCCGATCT",
                    "CBC": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 16},
                    "UMI": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 10},
                    "SLS": "TTTCTTATATGGG",  # Switch Leader Seq
                    "cDNA": RANDOM_SEGMENT_NAME,
                    "Poly_A": {HPR_SEGMENT_TYPE_NAME: ("A", 30)},
                    "3p_Adapter": "GTACTCTGCGTTGATACCACTGCTT",
                },
                "named_random_segments": ["CBC", "UMI", "cDNA"],
                "coding_region": "cDNA",
                "annotation_segments": {
                    "UMI": [(longbow.utils.constants.READ_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG),
                            (longbow.utils.constants.READ_RAW_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG)],
                    "CBC": [(longbow.utils.constants.READ_BARCODE_TAG, longbow.utils.constants.READ_BARCODE_POS_TAG),
                            (longbow.utils.constants.READ_RAW_BARCODE_TAG, longbow.utils.constants.READ_BARCODE_POS_TAG)],
                },
                "deprecated": False,
            },

            "bulk_10x5p_indexed": {
                "description": "bulk 10x 5' kit",
                "version": "3.0.0",
                "structure": [ "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter" ],
                "adapters": {
                    "5p_Adapter": "TCTACACGACGCTCTTCCGATCT",
                    "UMI": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 10},
                    "SLS": "TTTCTTATATGGG",  # Switch Leader Seq
                    "cDNA": RANDOM_SEGMENT_NAME,
                    "Poly_A": {HPR_SEGMENT_TYPE_NAME: ("A", 30)},
                    "sample_index": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 10},
                    "3p_Adapter": "CTCTGCGTTGATACCACTGCTT",
                },
                "named_random_segments": ["UMI", "cDNA", "sample_index"],
                "coding_region": "cDNA",
                "annotation_segments": {
                    "UMI": [(longbow.utils.constants.READ_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG),
                            (longbow.utils.constants.READ_RAW_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG)],
                    "sample_index": [(longbow.utils.constants.READ_DEMUX_TAG, longbow.utils.constants.READ_DEMUX_POS_TAG)],
                },
                "deprecated": False,
            },

            "bulk_10x5p_indexed_debug": {
                "description": "bulk 10x 5' kit",
                "version": "3.0.0",
                "structure": [ "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter",
                               "5p_Adapter'", "UMI", "sls", "cDNA", "Poly_T", "sample_index", "3p_Adapter'" ],
                "adapters": {
                    "5p_Adapter": "TCTACACGACGCTCTTCCGATCT",
                    "5p_Adapter'": "AGATCGGAAGAGCGTCGTGTAGA",
                    "UMI": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 12},
                    "SLS": "TTTCTTATATGGG",  # Switch Leader Seq
                    "sls": "CCCATATAAGAAA",  # Switch Leader Seq
                    "cDNA": RANDOM_SEGMENT_NAME,
                    "Poly_A": {HPR_SEGMENT_TYPE_NAME: ("A", 30)},
                    "Poly_T": {HPR_SEGMENT_TYPE_NAME: ("T", 30)},
                    "sample_index": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 10},
                    "3p_Adapter": "CTCTGCGTTGATACCACTGCTT",
                    "3p_Adapter'": "AAGCAGTGGTATCAACGCAGAG",
                },
                "named_random_segments": ["UMI", "cDNA", "sample_index"],
                "coding_region": "cDNA",
                "annotation_segments": {
                    "UMI": [(longbow.utils.constants.READ_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG),
                            (longbow.utils.constants.READ_RAW_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG)],
                    "sample_index": [(longbow.utils.constants.READ_DEMUX_TAG, longbow.utils.constants.READ_DEMUX_POS_TAG)],
                },
                "deprecated": False,
            },

            "bulk_teloprimeV2": {
                "description": "Lexogen TeloPrime V2 kit",
                "version": "3.0.0",
                "structure": [ "TPV2_adapter", "cDNA", "Poly_A", "idx", "rev_bind" ],
                "adapters": {
                    "TPV2_adapter": "CTACACGACGCTCTTCCGATCTTGGATTGATATGTAATACGACTCACTATAG",
                    "cDNA": RANDOM_SEGMENT_NAME,
                    "Poly_A": {HPR_SEGMENT_TYPE_NAME: ("A", 30)},
                    "idx": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 10},
                    "rev_bind": "CTCTGCGTTGATACCACTGCTT",
                },
                "named_random_segments": ["idx", "cDNA"],
                "coding_region": "cDNA",
                "annotation_segments": {
                    "idx": [(longbow.utils.constants.READ_INDEX_TAG, longbow.utils.constants.READ_BARCODE_POS_TAG)],
                },
                "deprecated": False,
            },

            # The slide-seq model is:
            #
            #                 |-----5p_Adapter---->        |--splitter------>               |------Poly_T---------------->                  |--------5p_Adapter----------|                     # noqa
            # AGCTTACTTGTGAAGACTACACGACGCTCTTCCGATCTNNNNNNNNTCTTCAGCGTTCCCGAGANNNNNNNNNNNNNVVTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVNNNNNNNNNNNNNNNNNCCCATGTACTCTGCGTTGATACCACTGCTTACTTGTAAGCTGTCTA... # noqa
            # |------A------->                      <------|                  <-----------|                                 <----cDNA-------|                              |-------B------>    # noqa
            #                                          V                           V
            #                                    Spatial Barcode 2         Spatial Barcode 1
            "spatial_slideseq": {
                "description": "Slide-seq protocol",
                "version": "3.0.0",
                "structure": [ "5p_Adapter", "SBC2", "SLS2", "SBC1", "UMI", "Poly_T", "cDNA", "3p_Adapter" ],
                "adapters": {
                    "5p_Adapter": "TCTACACGACGCTCTTCCGATCT",
                    "SBC2": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 8},
                    "SLS2": "TCTTCAGCGTTCCCGAGA",  # Switch Leader Seq
                    "SBC1": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 6},
                    # The UMI might be 7, rather than 9 elements long - not clear from the geneious file.
                    "UMI": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 9},
                    "Poly_T": {HPR_SEGMENT_TYPE_NAME: ("T", 30)},
                    "cDNA": RANDOM_SEGMENT_NAME,
                    "3p_Adapter": "CCCATGTACTCTGCGTTGATACCACTGCTT",
                },
                "named_random_segments": ["UMI", "SBC2", "SBC1", "cDNA"],
                "coding_region": "cDNA",
                "annotation_segments": {
                    "UMI": [(longbow.utils.constants.READ_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG),
                            (longbow.utils.constants.READ_RAW_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG)],
                    "SBC1": [(longbow.utils.constants.READ_SPATIAL_BARCODE1_TAG,
                            longbow.utils.constants.READ_SPATIAL_BARCODE1_POS_TAG)],
                    "SBC2": [(longbow.utils.constants.READ_SPATIAL_BARCODE2_TAG,
                            longbow.utils.constants.READ_SPATIAL_BARCODE2_POS_TAG)],
                },
                "deprecated": False,
            },
        },

        'multiplex': {
            "smrtbell_plate_3.0": {
                "description": "Multiplexing barcodes from SMRTbell barcoded adapter plate 3.0",
                "version": "3.0.0",
                "structure": [ "index" ],
                "adapters": {
                    "index": {
                        "bc2001": "ATCGTGCGACGAGTAT",
                        "bc2002": "TGCATGTCATGAGTAT",
                        "bc2003": "ACGAGTGCTCGAGTAT",
                        "bc2004": "TGCAGTGCTCGAGTAT",
                        "bc2005": "TGACTCGATCGAGTAT",
                        "bc2006": "CATGCGATCTGAGTAT",
                        "bc2007": "ACTAGCATCTGAGTAT",
                        "bc2008": "ACGCTAGTCTGAGTAT",
                        "bc2009": "CGATCGCACTGAGTAT",
                        "bc2010": "TACGTAGTATGAGTAT",
                        "bc2011": "CTGACAGTACGAGTAT",
                        "bc2012": "TCGTACTACTGAGTAT",
                        "bc2013": "CTGCGTAGACGAGTAT",
                        "bc2014": "ATACATGCACGAGTAT",
                        "bc2015": "CGACATAGATGAGTAT",
                        "bc2016": "ATCTGCACGTGAGTAT",
                        "bc2017": "CTATGATAGCGAGTAT",
                        "bc2018": "CGATCAGTGCGAGTAT",
                        "bc2019": "CGTCATAGTCGAGTAT",
                        "bc2020": "ACTATGCGTCGAGTAT",
                        "bc2021": "CGTACATGCTGAGTAT",
                        "bc2022": "TCATCGACGTGAGTAT",
                        "bc2023": "TCGCATGACTGAGTAT",
                        "bc2024": "CATGATCGACGAGTAT",
                        "bc2025": "ACGCACGTACGAGTAT",
                        "bc2026": "CAGTAGCGTCGAGTAT",
                        "bc2027": "TGACTGTAGCGAGTAT",
                        "bc2028": "ACTGCAGCACGAGTAT",
                        "bc2029": "TAGCAGTATCGAGTAT",
                        "bc2030": "CATACAGCATGAGTAT",
                        "bc2031": "ATAGCGTACTGAGTAT",
                        "bc2032": "ATAGACGAGTGAGTAT",
                        "bc2033": "CGACTCGTATGAGTAT",
                        "bc2034": "TACTAGTGACGAGTAT",
                        "bc2035": "CAGCTGACATGAGTAT",
                        "bc2036": "ACGTCGCTGCGAGTAT",
                        "bc2037": "CAGTATGAGCGAGTAT",
                        "bc2038": "TCACGACGACGAGTAT",
                        "bc2039": "CATGTATGTCGAGTAT",
                        "bc2040": "TGCTGCGACTGAGTAT",
                        "bc2041": "TATGATCACTGAGTAT",
                        "bc2042": "TCTGCACTGCGAGTAT",
                        "bc2043": "ACGATGACGTGAGTAT",
                        "bc2044": "CGATGATGCTGAGTAT",
                        "bc2045": "TACGACAGTCGAGTAT",
                        "bc2046": "ATAGCATGTCGAGTAT",
                        "bc2047": "CATAGTACTCGAGTAT",
                        "bc2048": "TGATGCTAGTGAGTAT",
                        "bc2049": "TAGTCTGCGTGAGTAT",
                        "bc2050": "CTCATCTATCGAGTAT",
                        "bc2051": "TGCATACTGCGAGTAT",
                        "bc2052": "CAGACTAGTCGAGTAT",
                        "bc2053": "ATCGTGATCTGAGTAT",
                        "bc2054": "CTGCGATCACGAGTAT",
                        "bc2055": "CTCAGCATACGAGTAT",
                        "bc2056": "TCGCAGCGTCGAGTAT",
                        "bc2057": "TAGCACGCATGAGTAT",
                        "bc2058": "TACTGACGCTGAGTAT",
                        "bc2059": "ATCTGACTATGAGTAT",
                        "bc2060": "ATACGAGCTCGAGTAT",
                        "bc2061": "CGAGCACGCTGAGTAT",
                        "bc2062": "TCTGCGTATCGAGTAT",
                        "bc2063": "TCTGCATCATGAGTAT",
                        "bc2064": "TGCGTGATGCGAGTAT",
                        "bc2065": "TGAGCTATGCGAGTAT",
                        "bc2066": "CTGTCGTAGTGAGTAT",
                        "bc2067": "ATCGATGCATGAGTAT",
                        "bc2068": "ACTACGTGATGAGTAT",
                        "bc2069": "TCTATGACATGAGTAT",
                        "bc2070": "TACTGCTCACGAGTAT",
                        "bc2071": "CGAGTCTAGCGAGTAT",
                        "bc2072": "TATCAGTAGTGAGTAT",
                        "bc2073": "ATCACTAGTCGAGTAT",
                        "bc2074": "TATCACGACTGAGTAT",
                        "bc2075": "CTCGTCAGATGAGTAT",
                        "bc2076": "CAGCAGTGACGAGTAT",
                        "bc2077": "TGCGACGTGCGAGTAT",
                        "bc2078": "CTCACTGAGTGAGTAT",
                        "bc2079": "CACTGAGCGTGAGTAT",
                        "bc2080": "CAGCGTCTACGAGTAT",
                        "bc2081": "CTACTATGTCGAGTAT",
                        "bc2082": "ATGTACAGACGAGTAT",
                        "bc2083": "ACTCATCAGTGAGTAT",
                        "bc2084": "CTGAGCACTCGAGTAT",
                        "bc2085": "ATCATCTACTGAGTAT",
                        "bc2086": "TACATGCGATGAGTAT",
                        "bc2087": "TCGCTGTCACGAGTAT",
                        "bc2088": "ACGCTCATGCGAGTAT",
                        "bc2089": "TACTAGCAGCGAGTAT",
                        "bc2090": "CGTAGCAGATGAGTAT",
                        "bc2091": "CGTGCTCGTCGAGTAT",
                        "bc2092": "ACAGCTGTACGAGTAT",
                        "bc2093": "TCGATGCTACGAGTAT",
                        "bc2094": "TAGATACAGCGAGTAT",
                        "bc2095": "CTACTCATACGAGTAT",
                        "bc2096": "ATGTACTAGTGAGTAT"
                    }
                },
                "deprecated": False,
            },
        }
    }