import re
import json
import logging

import click_log

import numpy as np
import pandas as pd

from pomegranate import *
from pomegranate.callbacks import History, ModelCheckpoint


logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger(__name__)
click_log.basic_config(logger)


class LibraryModel:
    """Model describing a given library preparation.
    The model can annotate the known sections of a read from the library it describes."""

    # Define constants for all our default probabilities here:
    RANDOM_BASE_PROB = 0.5
    PER_BASE_MATCH_PROB = 0.94
    PER_BASE_SNP_PROB = 0.02

    MATCH_MATCH_PROB = 0.90
    MATCH_INDEL_PROB = 0.05
    MATCH_TRAIL_INSERT_PROB = 0.90

    INDEL_CONTINUATION_PROB = 0.7
    INDEL_SWITCH_PROB = 0.3

    START_AND_END_RANDOM_PROB = 0.5
    RAND_RAND_PROB = 0.5
    RAND_INS_CONTINUATION_PROB = 0.8
    RAND_INS_TO_DEL_PROB = 0.1
    RAND_INS_END_PROB = 0.1

    SUDDEN_END_PROB = 0.01
    MATCH_END_PROB = 0.1

    def __init__(self,
                 name,
                 array_element_structure,
                 adapters,
                 direct_connections,
                 start_element_names,
                 end_element_names,
                 do_build=True):

        self.name = name

        self.array_element_structure = array_element_structure
        self.adapter_dict = adapters
        self.direct_connections_dict = direct_connections

        self.start_element_names = start_element_names
        self.end_element_names = end_element_names

        self.hmm = None
        self.key_adapters = self._create_key_adapter_order()
        self.key_adapter_set = set(self.key_adapters)

        if do_build:
            self.build()

    def annotate(self, seq):
        """Annotate the given segment using this model."""
        logp, path = self.hmm.viterbi(seq)

        ppath = []
        for p, (idx, state) in enumerate(path[1:-1]):
            if (
                    "start" not in state.name
                    and ":RD" not in state.name
                    and ":D" not in state.name
            ):
                ppath.append(f'{re.split(":", state.name)[0]}')

        return logp, ppath

    def validate_segment_order(self, ordered_segment_names, allow_missing_first_adapter=True):
        """Validate the order of the given segments against the expected order in this model.

        Returns: (True|False, # key adapters found, first key adapter index)"""

        # Iterate through our given segment names and check if they occur in order:
        num_key_adapters_found = 0
        key_adapter_indx = 0
        first_key_adapter_index = 0
        found_first_key = False
        for n in ordered_segment_names:
            # Ignore all segment names that do not characterize our library:
            if n in self.key_adapter_set:

                # If this is our first segment, we should allow for the possibility that our array begins
                # somewhere after the first element.  We must find the starting point:
                if not found_first_key:
                    while n != self.key_adapters[key_adapter_indx]:
                        key_adapter_indx += 1
                    first_key_adapter_index = key_adapter_indx
                    found_first_key = True

                    # TODO: This can be eliminated for newer datasets, but is here because we're still testing with
                    #       older data that does not have the first MAS-seq overhang.  The model has already been
                    #       updated to handle these segments.
                    if allow_missing_first_adapter and (key_adapter_indx == 1) and \
                            (ordered_segment_names[0] == "10x_Adapter"):
                        num_key_adapters_found += 1
                        first_key_adapter_index = 0

                # Check our key segments here:
                # NOTE: it's possible that we can start matching in the middle of an array, but at the start of a read,
                #       and so we have to do a bounds check here:
                if key_adapter_indx < len(self.key_adapters) and (n == self.key_adapters[key_adapter_indx]):
                    key_adapter_indx += 1
                    num_key_adapters_found += 1
                else:
                    # This read does not conform to the model!
                    return False, num_key_adapters_found, first_key_adapter_index

        # If we've made it here and we have seen at least 1 key adapter, then we have a valid array:
        is_valid = True if num_key_adapters_found > 0 else False
        return is_valid, num_key_adapters_found, first_key_adapter_index

    def extract_key_segment_names(self, segment_names):
        """Return a list of key segment names from the given list of segment_names."""
        return [n for n in segment_names if n in self.key_adapter_set]

    def build(self):
        """Build the HMM underlying this model given our segment information."""
        self.hmm = LibraryModel._make_random_repeat_model()
        for k, v in self.adapter_dict.items():
            self.hmm.add_model(LibraryModel._make_global_alignment_model(v, k))

        self.hmm.bake(merge="None")

        # dictionary of model starting states, random start, and random end
        starts = {}
        rda = None
        rdb = None
        for s in self.hmm.states:
            if "-start" in s.name and "random" not in s.name:
                starts[re.sub("-start", "", s.name)] = s
            elif "random:RDA" in s.name:
                rda = s
            elif "random:RDB" in s.name:
                rdb = s

        # link array element start to hmm start nodes:
        for sname in starts:
            if sname in self.start_element_names:
                self.hmm.add_transition(self.hmm.start, starts[sname], 1.0/len(self.start_element_names))

        # link array element ends to start nodes:
        for sname in starts:
            if sname in self.end_element_names:
                self.hmm.add_transition(rda, starts[sname], 1.0/len(self.end_element_names))

        # link up ending states according to our direct connections dictionary
        for s in self.hmm.states:
            m = re.match(r"^(\w+):([MID])(\d+)", s.name)
            if m is not None and int(m.group(3)) == len(self.adapter_dict[m.group(1)]):
                sname = m.group(1)

                if sname in self.direct_connections_dict:
                    for dcname in self.direct_connections_dict[sname]:
                        self.hmm.add_transition(
                            s, starts[dcname], 1.0 / len(self.direct_connections_dict[sname])
                        )
                else:
                    # Verify this probability is the notional equivalent:
                    # self.hmm.add_transition(s, rdb, 0.5)
                    self.hmm.add_transition(s, rdb, LibraryModel.RAND_RAND_PROB)

        # link up all adapters to model end state
        for s in self.hmm.states:
            m = re.match(r"^(\w+):([MID])(\d+)", s.name)
            if m is not None and int(m.group(3)) == len(self.adapter_dict[m.group(1)]):
                self.hmm.add_transition(s, self.hmm.end, LibraryModel.SUDDEN_END_PROB)

        self.hmm.bake()

    def _create_key_adapter_order(self):
        """Setup an ordered list of key segments that characterize the correct array element order."""

        # TODO: Generalize this for all library types / segment names!
        # Assumption: self.array_element_structure contains the array elements in our library in the order in which
        #             they appear in the data.
        # Heuristic: The segments that characterize the array elements themselves all have single-character names, so we
        #            filter out all segments from self.array_element_structure with names longer than 1 char.  We then
        #            use these in order to characterize the reads.

        ordered_key_adapters = [s for array in self.array_element_structure for s in array if len(s) == 1]
        return ordered_key_adapters

    @staticmethod
    def _make_global_alignment_model(target, name=None):
        model = HiddenMarkovModel(name=name)
        s = {}

        # add states
        i0 = State(
            DiscreteDistribution({
                    "A": LibraryModel.RANDOM_BASE_PROB,
                    "C": LibraryModel.RANDOM_BASE_PROB,
                    "G": LibraryModel.RANDOM_BASE_PROB,
                    "T": LibraryModel.RANDOM_BASE_PROB
                }),
            name=f"{name}:I0",
        )

        model.add_state(i0)

        s[i0.name] = i0

        for c in range(len(target)):
            dc = State(None, name=f"{name}:D{c + 1}")

            mc = State(
                DiscreteDistribution(
                    {
                        "A": LibraryModel.PER_BASE_MATCH_PROB if target[c] == "A" else LibraryModel.PER_BASE_SNP_PROB,
                        "C": LibraryModel.PER_BASE_MATCH_PROB if target[c] == "C" else LibraryModel.PER_BASE_SNP_PROB,
                        "G": LibraryModel.PER_BASE_MATCH_PROB if target[c] == "G" else LibraryModel.PER_BASE_SNP_PROB,
                        "T": LibraryModel.PER_BASE_MATCH_PROB if target[c] == "T" else LibraryModel.PER_BASE_SNP_PROB,
                    }
                ),
                name=f"{name}:M{c + 1}",
            )

            ic = State(
                DiscreteDistribution({
                    "A": LibraryModel.RANDOM_BASE_PROB,
                    "C": LibraryModel.RANDOM_BASE_PROB,
                    "G": LibraryModel.RANDOM_BASE_PROB,
                    "T": LibraryModel.RANDOM_BASE_PROB
                }),
                name=f"{name}:I{c + 1}",
            )

            model.add_states([mc, ic, dc])

            s[dc.name] = dc
            s[mc.name] = mc
            s[ic.name] = ic

        # add transitions
        model.add_transition(model.start, s[f"{name}:I0"], LibraryModel.MATCH_INDEL_PROB)
        model.add_transition(model.start, s[f"{name}:D1"], LibraryModel.MATCH_INDEL_PROB)
        model.add_transition(model.start, s[f"{name}:M1"], LibraryModel.MATCH_MATCH_PROB)

        model.add_transition(s[f"{name}:I0"], s[f"{name}:I0"], LibraryModel.INDEL_CONTINUATION_PROB)
        model.add_transition(s[f"{name}:I0"], s[f"{name}:D1"], LibraryModel.INDEL_SWITCH_PROB/2)
        model.add_transition(s[f"{name}:I0"], s[f"{name}:M1"], LibraryModel.INDEL_SWITCH_PROB/2)

        for c in range(1, len(target)):
            model.add_transition(s[f"{name}:D{c}"], s[f"{name}:D{c + 1}"], LibraryModel.INDEL_SWITCH_PROB/2)
            model.add_transition(s[f"{name}:D{c}"], s[f"{name}:I{c}"], LibraryModel.INDEL_CONTINUATION_PROB)
            model.add_transition(s[f"{name}:D{c}"], s[f"{name}:M{c + 1}"], LibraryModel.INDEL_SWITCH_PROB/2)

            model.add_transition(s[f"{name}:I{c}"], s[f"{name}:D{c + 1}"], LibraryModel.INDEL_SWITCH_PROB/2)
            model.add_transition(s[f"{name}:I{c}"], s[f"{name}:I{c}"], LibraryModel.INDEL_SWITCH_PROB/2)
            model.add_transition(s[f"{name}:I{c}"], s[f"{name}:M{c + 1}"], LibraryModel.INDEL_CONTINUATION_PROB)

            model.add_transition(s[f"{name}:M{c}"], s[f"{name}:D{c + 1}"], LibraryModel.MATCH_INDEL_PROB)
            model.add_transition(s[f"{name}:M{c}"], s[f"{name}:I{c}"], LibraryModel.MATCH_INDEL_PROB)
            model.add_transition(s[f"{name}:M{c}"], s[f"{name}:M{c + 1}"], LibraryModel.MATCH_MATCH_PROB)

        model.add_transition(
            s[f"{name}:D{len(target)}"], s[f"{name}:I{len(target)}"], LibraryModel.INDEL_CONTINUATION_PROB
        )
        model.add_transition(s[f"{name}:D{len(target)}"], model.end, LibraryModel.INDEL_SWITCH_PROB)

        model.add_transition(s[f"{name}:I{len(target)}"], s[f"{name}:I{len(target)}"], LibraryModel.INDEL_SWITCH_PROB/2)
        model.add_transition(
            s[f"{name}:I{len(target)}"],
            model.end,
            LibraryModel.INDEL_CONTINUATION_PROB + (LibraryModel.INDEL_SWITCH_PROB/2)
        )

        model.add_transition(
            s[f"{name}:M{len(target)}"], s[f"{name}:I{len(target)}"], LibraryModel.MATCH_TRAIL_INSERT_PROB
        )
        model.add_transition(s[f"{name}:M{len(target)}"], model.end, LibraryModel.MATCH_END_PROB)

        model.bake(merge="None")

        return model

    @staticmethod
    def _make_random_repeat_model(name="random"):
        model = HiddenMarkovModel(name=name)

        # add states
        ri = State(
            DiscreteDistribution({
                "A": LibraryModel.RANDOM_BASE_PROB,
                "C": LibraryModel.RANDOM_BASE_PROB,
                "G": LibraryModel.RANDOM_BASE_PROB,
                "T": LibraryModel.RANDOM_BASE_PROB
            }),
            name=f"{name}:RI",
        )
        rda = State(None, name=f"{name}:RDA")
        rdb = State(None, name=f"{name}:RDB")

        model.add_states([ri, rda, rdb])

        # add transitions
        model.add_transition(model.start, rda, LibraryModel.START_AND_END_RANDOM_PROB)
        model.add_transition(model.start, ri, LibraryModel.START_AND_END_RANDOM_PROB)

        model.add_transition(ri, ri, LibraryModel.RAND_INS_CONTINUATION_PROB)
        model.add_transition(ri, rda, LibraryModel.RAND_INS_TO_DEL_PROB)
        model.add_transition(ri, model.end, LibraryModel.RAND_INS_END_PROB)

        model.add_transition(rdb, ri, LibraryModel.RAND_RAND_PROB)
        model.add_transition(rdb, model.end, LibraryModel.START_AND_END_RANDOM_PROB)

        model.bake(merge="None")

        return model

    # TODO: FINISH THIS!
    def to_json(self, outfile=None, indent=4):
        """Serialize this model to a json object and return that json object.
        If outfile is not none, will write the json object to the given file path."""

        model_data = {
            "name": self.name,
            "array_element_structure": self.array_element_structure,
            "adapters": self.adapter_dict,
            "direct_connections": {k: list(v) for k, v in self.direct_connections_dict.items()},
            "start_element_names": list(self.start_element_names),
            "end_element_names": list(self.end_element_names)
        }

        if outfile:
            with open(outfile, 'w') as f:
                json.dump(model_data, f, indent=indent)

        return json.dumps(model_data, indent=indent)

    # TODO: FINISH THIS!
    @staticmethod
    def from_json_file(json_file):
        """Create a LibraryModel instance from the given json data.
        This method will open the file at the given location and use the data in that file to create a LibraryModel."""

        try:
            with open(json_file) as f:
                json_data = json.load(f)
        except FileNotFoundError:
            logger.error(f"File does not exist: {json_file}")
            sys.exit(1)

        m = LibraryModel(
            name=json_data["name"],
            array_element_structure=tuple(tuple(v) for v in json_data["array_element_structure"]),
            adapters=json_data["adapters"],
            direct_connections={k: set(v) for k, v in json_data["direct_connections"].items()},
            start_element_names=set(json_data["start_element_names"]),
            end_element_names=set(json_data["end_element_names"]),
        )

        print(m.name)
        print(m.array_element_structure)
        print(m.adapter_dict)
        print(m.direct_connections_dict)
        print(m.start_element_names)
        print(m.end_element_names)

        return m

    @staticmethod
    def build_and_return_mas_seq_model():
        """Create and return the model for the standard 15 element MAS-seq array."""
        return LibraryModel(
            name="mas15",
            array_element_structure=(
                # NOTE: the first element doesn't currently have the "A" adapter in this version of the library.
                ("A", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("B", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("C", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("D", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("E", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("F", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("G", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("H", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("I", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("J", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("K", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("L", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("M", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("N", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                # The last element doesn't currently have the "P" adapter in this version of the library:
                ("O", "10x_Adapter", "random", "Poly_A", "3p_Adapter", "P"),
            ),
            adapters={
                "10x_Adapter": "TCTACACGACGCTCTTCCGATCT",
                "Poly_A": "A" * 30,
                "3p_Adapter": "GTACTCTGCGTTGATACCACTGCTT",
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
            direct_connections={
                "Poly_A": {"3p_Adapter"},
                "3p_Adapter": {
                    "A",
                    "B",
                    "C",
                    "D",
                    "E",
                    "F",
                    "G",
                    "H",
                    "I",
                    "J",
                    "K",
                    "L",
                    "M",
                    "N",
                    "O",
                    "P",
                },
                "A": {"10x_Adapter"},
                "B": {"10x_Adapter"},
                "C": {"10x_Adapter"},
                "D": {"10x_Adapter"},
                "E": {"10x_Adapter"},
                "F": {"10x_Adapter"},
                "G": {"10x_Adapter"},
                "H": {"10x_Adapter"},
                "I": {"10x_Adapter"},
                "J": {"10x_Adapter"},
                "K": {"10x_Adapter"},
                "L": {"10x_Adapter"},
                "M": {"10x_Adapter"},
                "N": {"10x_Adapter"},
                "O": {"10x_Adapter"},
                "P": {"10x_Adapter"},
            },
            start_element_names={"A", "10x_Adapter"},
            end_element_names={"Poly_A", "P"},
        )

    @staticmethod
    def build_and_return_mas_seq_10_model():
        """Create and return the model for the 10 element MAS-seq array."""
        return LibraryModel(
            name="mas10",
            array_element_structure=(
                ("Q", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("C", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("M", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("I", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("O", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("J", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("B", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("D", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                ("K", "10x_Adapter", "random", "Poly_A", "3p_Adapter"),
                # The last element may not currently have the "R" adapter in this version of the library:
                ("H", "10x_Adapter", "random", "Poly_A", "3p_Adapter", "R"),
            ),
            adapters={
                "10x_Adapter": "TCTACACGACGCTCTTCCGATCT",
                "Poly_A": "A" * 30,
                "3p_Adapter": "GTACTCTGCGTTGATACCACTGCTT",
                "B": "ACTTGTAAGCTGTCTA",
                "C": "ACTCTGTCAGGTCCGA",
                "D": "ACCTCCTCCTCCAGAA",
                "H": "ATGTTGAATCCTAGCG",
                "I": "AGTGCGTTGCGAATTG",
                "J": "AATTGCGTAGTTGGCC",
                "K": "ACACTTGGTCGCAATC",
                "M": "ACCTAGATCAGAGCCT",
                "O": "AAGTCACCGGCACCTT",
                "Q": "AAGCACCATAATGTGT",
                "R": "AACCGGACACACTTAG",
            },
            direct_connections={
                "Poly_A": {"3p_Adapter"},
                "3p_Adapter": {
                    "Q",
                    "C",
                    "M",
                    "I",
                    "O",
                    "J",
                    "B",
                    "D",
                    "K",
                    "H",
                    "R",
                },
                "B": {"10x_Adapter"},
                "C": {"10x_Adapter"},
                "D": {"10x_Adapter"},
                "H": {"10x_Adapter"},
                "I": {"10x_Adapter"},
                "J": {"10x_Adapter"},
                "K": {"10x_Adapter"},
                "M": {"10x_Adapter"},
                "O": {"10x_Adapter"},
                "Q": {"10x_Adapter"},
                "R": {"10x_Adapter"},
            },
            start_element_names={"Q", "10x_Adapter"},
            end_element_names={"Poly_A", "R"},
        )


# IUPAC RC's from: http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
# and https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
RC_BASE_MAP = {
    "N": "N",
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "Y": "R",
    "R": "Y",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "V": "B",
    "D": "H",
    "H": "D",
    "n": "n",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
    "y": "r",
    "r": "y",
    "s": "s",
    "w": "w",
    "k": "m",
    "m": "k",
    "b": "v",
    "v": "b",
    "d": "h",
    "h": "d",
}


def reverse_complement(base_string):
    """
    Reverse complements the given base_string.
    :param base_string: String of bases to be reverse-complemented.
    :return: The reverse complement of the given base string.
    """

    return "".join(map(lambda b: RC_BASE_MAP[b], base_string[::-1]))
