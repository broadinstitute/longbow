import sys
import re
import importlib.resources
import json
import logging

import click_log

from pomegranate import *

from .constants import (
    RANDOM_SEGMENT_NAME,
    RANDOM_SILENT_STATE_A,
    RANDOM_SILENT_STATE_B,
    RANDOM_BASE_STATE,
    BAKE_MERGE_STRATEGY,
)

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger(__name__)
click_log.basic_config(logger)

starts_with_number_re = re.compile(r"^\d")


def load_preconfigured_models():
    pre_configured_models = {"array": {}, "cdna": {}}

    with importlib.resources.path("longbow", "preconfigured_models") as model_dir:
        for json_file in (model_dir / "array").glob("*json"):
            with json_file.open() as fh:
                m = json.load(fh)
                pre_configured_models["array"][m["name"]] = m

        for json_file in (model_dir / "cdna").glob("*json"):
            with json_file.open() as fh:
                m = json.load(fh)
                pre_configured_models["cdna"][m["name"]] = m

    return pre_configured_models


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

    pre_configured_models = load_preconfigured_models()

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
        model.add_transition(s[f"{name}:I{expected_length}"], s[f"{name}:I{expected_length}"],
                             ModelBuilder.INDEL_SWITCH_PROB / 2)
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
        logger.debug("Making Model: RANDOM REPEAT")
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
