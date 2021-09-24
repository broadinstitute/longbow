import sys
import re
import json
import logging
import queue

import networkx

import click_log

import numpy as np
import pandas as pd

from pomegranate import *
from pomegranate.callbacks import History, ModelCheckpoint


logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger(__name__)
click_log.basic_config(logger)

RANDOM_SEGMENT_NAME = "random"
FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME = "FixedLengthRandomBases"


class LibraryModel:
    """Model describing a given library preparation.
    The model can annotate the known sections of a read from the library it describes."""

    # Define constants for all our default probabilities here:
    RANDOM_BASE_PROB = 0.5
    PER_BASE_MATCH_PROB = 0.94
    PER_BASE_MISMATCH_PROB = 0.02

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
                 description,
                 version,
                 array_element_structure,
                 adapters,
                 direct_connections,
                 start_element_names,
                 end_element_names,
                 named_random_segments,
                 do_build=True):

        self.name = name
        self.description = description
        self.version = version

        self.array_element_structure = array_element_structure
        self.adapter_dict = adapters
        self.direct_connections_dict = direct_connections

        self.start_element_names = start_element_names
        self.end_element_names = end_element_names

        self.named_random_segments = named_random_segments

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

        # Except!  If we have an array singleton and it's the terminal overhang adapter,
        # then we do NOT have a valid array:
        if is_valid and (num_key_adapters_found == 1) and \
                (self.key_adapters[first_key_adapter_index] == self.array_element_structure[-1][-1]):
            is_valid = False

        return is_valid, num_key_adapters_found, first_key_adapter_index

    def extract_key_segment_names(self, segment_names):
        """Return a list of key segment names from the given list of segment_names."""
        return [n for n in segment_names if n in self.key_adapter_set]

    # def build(self):
    #     """Build the HMM underlying this model given our segment information."""
    #     self.hmm = LibraryModel._make_random_repeat_model()
    #     for k, v in self.adapter_dict.items():
    #         self.hmm.add_model(LibraryModel._make_global_alignment_model(v, k))
    #
    #     self.hmm.bake(merge="None")
    #
    #     # dictionary of model starting states, random start, and random end
    #     starts = {}
    #     rda = None
    #     rdb = None
    #     for s in self.hmm.states:
    #         if "-start" in s.name and RANDOM_SEGMENT_NAME not in s.name:
    #             starts[re.sub("-start", "", s.name)] = s
    #         elif s.name.startswith(f"{RANDOM_SEGMENT_NAME}:RDA"):
    #             rda = s
    #         elif s.name.startswith(f"{RANDOM_SEGMENT_NAME}:RDB"):
    #             rdb = s
    #
    #     # link silent start state to all of our start nodes
    #     # this is the start of the array.
    #     for sname in starts:
    #         if sname in self.start_element_names:
    #             self.hmm.add_transition(self.hmm.start, starts[sname], 1.0 / len(self.start_element_names))
    #
    #     # link random element model silent state A to all of our starts
    #     # This is the transition from `random` to one of our array elements.
    #     for sname in starts:
    #         if sname in self.end_element_names:
    #             self.hmm.add_transition(rda, starts[sname], 1.0 / len(self.end_element_names))
    #
    #     # link up adapter final states according to our direct connections dictionary
    #     for s in self.hmm.states:
    #         m = re.match(r"^(\w+):([MID])(\d+)", s.name)
    #
    #         # If we're in the last state of each adapter:
    #         if m is not None and int(m.group(3)) == len(self.adapter_dict[m.group(1)]):
    #             sname = m.group(1)
    #
    #             # If we have direct connections from this state to some other states then we add
    #             # them into the model:
    #             if sname in self.direct_connections_dict:
    #                 for dcname in self.direct_connections_dict[sname]:
    #                     self.hmm.add_transition(
    #                         s, starts[dcname], 1.0 / len(self.direct_connections_dict[sname])
    #                     )
    #             # If we have no direct connections from this state, we add in a transition from this
    #             # state to the silent random B state:
    #             else:
    #                 # Verify this probability is the notional equivalent:
    #                 self.hmm.add_transition(s, rdb, LibraryModel.RAND_RAND_PROB)
    #
    #     # link up all adapters to model end state
    #     for s in self.hmm.states:
    #         m = re.match(r"^(\w+):([MID])(\d+)", s.name)
    #
    #         # Get the last state in each adapter:
    #         if m is not None and int(m.group(3)) == len(self.adapter_dict[m.group(1)]):
    #             self.hmm.add_transition(s, self.hmm.end, LibraryModel.SUDDEN_END_PROB)
    #
    #     self.hmm.bake()
    #
    # def build2(self):
    #     """Build the HMM underlying this model given our segment information."""
    #
    #     self.hmm = None
    #
    #     # Start with a random model, which will represent erroneous transitions between
    #     # self.hmm = LibraryModel._make_random_repeat_model(RANDOM_SEGMENT_NAME)
    #
    #     for k, v in self.adapter_dict.items():
    #         adapter_hmm = None
    #         if type(v) is str:
    #             if v in self.named_random_segments:
    #                 # We have a named random segment.
    #                 # We need to make a named random model:
    #                 adapter_hmm = LibraryModel._make_random_repeat_model(name=k)
    #             else:
    #                 # This must be a normal string of bases to model:
    #                 adapter_hmm = LibraryModel._make_global_alignment_model(v, k)
    #         elif type(v) is dict:
    #             segment_type = list(v.keys())[0]
    #
    #             # NOTE: Must add in new model types here when we add them:
    #             if segment_type == FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME:
    #                 # Fixed length random segments for barcodes (and similar):
    #                 adapter_hmm = LibraryModel._make_fixed_length_random_segment(k, list(v.values())[0])
    #             else:
    #                 logger.critical(f"Unknown special model type: {segment_type}")
    #                 raise RuntimeError(f"Unknown special model type: {segment_type}")
    #
    #         # Add this adapter's model to the HMM:
    #         if self.hmm is None:
    #             self.hmm = adapter_hmm
    #         else:
    #             self.hmm.add_model(adapter_hmm)
    #
    #     self.hmm.bake(merge="None")
    #
    #     # dictionary of model starting states, random start, and random end
    #     starts = {}
    #     rda = None
    #     rdb = None
    #     for s in self.hmm.states:
    #         state_name = LibraryModel._get_state_name(s)
    #         if s.name.endswith("-start") and state_name is not RANDOM_SEGMENT_NAME:
    #             starts[re.sub("-start", "", s.name)] = s
    #         elif s.name.startswith(f"{RANDOM_SEGMENT_NAME}:RDA"):
    #             rda = s
    #         elif s.name.startswith(f"{RANDOM_SEGMENT_NAME}:RDB"):
    #             rdb = s
    #
    #     # link silent start state to all of our start nodes
    #     # this is the start of the array.
    #     for sname in starts:
    #         if sname in self.start_element_names:
    #
    #             # Only add a transition of the states are not the same:
    #             if self.hmm.start == starts[sname]:
    #                 logger.warning(f"Ignoring attempt to add a simple cycle to graph: {starts[sname].name}->{starts[sname].name}")
    #                 continue
    #
    #             self.hmm.add_transition(self.hmm.start, starts[sname], 1.0/len(self.start_element_names))
    #
    #     # # link random element model silent state A to all of our starts
    #     # # This is the transition from `random` to one of our array elements.
    #     # for sname in starts:
    #     #     if sname in self.end_element_names:
    #     #         self.hmm.add_transition(rda, starts[sname], 1.0/len(self.end_element_names))
    #
    #     # link up adapter final states according to our direct connections dictionary
    #     state_type_regex = re.compile(r"^(\w+):([MID])(\d+)")
    #     for s in self.hmm.states:
    #         m = state_type_regex.match(s.name)
    #
    #         # If we're in the last state of each adapter:
    #         if m is not None and int(m.group(3)) == len(self.adapter_dict[m.group(1)]):
    #             sname = m.group(1)
    #
    #             # If we have direct connections from this state to some other states then we add
    #             # them into the model:
    #             if sname in self.direct_connections_dict:
    #                 for dcname in self.direct_connections_dict[sname]:
    #
    #                     # Only add a transition of the states are not the same:
    #                     if s == starts[dcname]:
    #                         logger.warning(
    #                             f"Ignoring attempt to add a simple cycle to graph: {starts[dcname].name}->{starts[dcname].name}")
    #                         continue
    #
    #                     self.hmm.add_transition(
    #                         s, starts[dcname], 1.0 / len(self.direct_connections_dict[sname])
    #                     )
    #             # # If we have no direct connections from this state, we add in a transition from this
    #             # # state to the silent random B state:
    #             # else:
    #             #     # Verify this probability is the notional equivalent:
    #             #     self.hmm.add_transition(s, rdb, LibraryModel.RAND_RAND_PROB)
    #
    #     # link up all adapters to model end state
    #     for s in self.hmm.states:
    #         m = state_type_regex.match(s.name)
    #
    #         # Get the last state in each adapter:
    #         if m is not None and int(m.group(3)) == len(self.adapter_dict[m.group(1)]):
    #             self.hmm.add_transition(s, self.hmm.end, LibraryModel.SUDDEN_END_PROB)
    #
    #     self.hmm.bake()

    def get_all_node_names(self):
        """Return a list containing the names of all nodes in this model."""
        all_node_names = {n for element in self.array_element_structure for n in element}
        for n in self.named_random_segments:
            all_node_names.add(n)
        for n in self.start_element_names:
            all_node_names.add(n)
        for n in self.end_element_names:
            all_node_names.add(n)
        for n in self.direct_connections_dict:
            all_node_names.add(n)
        for k in self.direct_connections_dict:
            for node_list in self.direct_connections_dict[k]:
                for n in node_list:
                    all_node_names.add(n)

        return sorted(list(all_node_names))

    def build(self):
        """Build the HMM underlying this model given our segment information."""

        # Validate our model here so we can go through and just worry about creating it:
        self.validate_model()

        # Create the centralized random repeat model that allows for errors between the known markers:
        self.hmm = LibraryModel._make_random_repeat_model(name=RANDOM_SEGMENT_NAME)

        # Create models for individual adapter sequences,
        # taking special care with named random sections:
        for k, v in self.adapter_dict.items():
            adapter_hmm = None
            if type(v) is str:
                if v in self.named_random_segments:
                    # We have a named random segment.
                    # We need to make a named random model:
                    adapter_hmm = LibraryModel._make_random_repeat_model(name=k)
                else:
                    # This must be a normal string of bases to model:
                    adapter_hmm = LibraryModel._make_global_alignment_model(v, k)
            elif type(v) is dict:
                segment_type = list(v.keys())[0]

                # NOTE: Must add in new model types here when we add them:
                if segment_type == FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME:
                    # Fixed length random segments for barcodes (and similar):
                    adapter_hmm = LibraryModel._make_fixed_length_random_segment(k, list(v.values())[0])
                else:
                    logger.critical(f"Unknown special model type: {segment_type}")
                    raise RuntimeError(f"Unknown special model type: {segment_type}")
            self.hmm.add_model(adapter_hmm)

        # Finalize all models that we have for now:
        self.hmm.bake(merge="None")

        # dictionary of model starting states, random start, and random end
        segment_starts = {}
        rda = None
        rdb = None
        for s in self.hmm.states:
            if "-start" in s.name and not s.name.startswith(RANDOM_SEGMENT_NAME):
                segment_starts[re.sub("-start", "", s.name)] = s
            elif s.name.startswith(f"{RANDOM_SEGMENT_NAME}:RDA"):
                rda = s
            elif s.name.startswith(f"{RANDOM_SEGMENT_NAME}:RDB"):
                rdb = s

        # link silent start state to all of our start nodes.
        # this is the start of the array.
        # NOTE: we cannot start in a named random segment.
        for sname in segment_starts:
            if sname in self.start_element_names:
                self.hmm.add_transition(self.hmm.start, segment_starts[sname], 1.0 / len(self.start_element_names))

        # link random element model silent state A to all of our starts
        # This is the transition from `random` to one of our array elements.
        # NOTE: we cannot end in a named random segment.
        for sname in segment_starts:
            if sname in self.end_element_names:
                self.hmm.add_transition(rda, segment_starts[sname], 1.0 / len(self.end_element_names))

        # # link random element model silent state A to all of our starts
        # # This is the transition from `random` to one of our array elements.
        # for sname in starts:
        #     num_standard_states = 0
        #     if sname not in self.named_random_segments:
        #         num_standard_states += 1
        # for sname in starts:
        #     self.hmm.add_transition(self.hmm.start, starts[sname], 1.0 / num_standard_states)

        # link up adapter final states according to our direct connections dictionary
        for s in self.hmm.states:
            m = re.match(r"^(\w+):([MID])(\d+)", s.name)

            # If we're in the last state of each adapter:
            if m is not None and int(m.group(3)) == len(self.adapter_dict[m.group(1)]):
                sname = m.group(1)

                # If we have direct connections from this state to some other states then we add
                # them into the model:
                if sname in self.direct_connections_dict:
                    for dcname in self.direct_connections_dict[sname]:
                        self.hmm.add_transition(
                            s, segment_starts[dcname], 1.0 / len(self.direct_connections_dict[sname])
                        )

                    # Also add in a connection to the silent random B state if this is a non-random section and it does
                    # not connect to a random section:
                    if sname not in self.named_random_segments:
                        is_connected_to_random = False
                        for n in self.direct_connections_dict[sname]:
                            if n in self.named_random_segments:
                                is_connected_to_random = True
                                break
                        if not is_connected_to_random:
                            # Verify this probability is the notional equivalent:
                            self.hmm.add_transition(s, rdb, LibraryModel.RAND_RAND_PROB)

                # If we have no direct connections from this state, we add in a transition from this
                # state to the silent random B state:
                else:
                    # Verify this probability is the notional equivalent:
                    self.hmm.add_transition(s, rdb, LibraryModel.RAND_RAND_PROB)

        # link up all adapters to model end state
        for s in self.hmm.states:
            m = re.match(r"^(\w+):([MID])(\d+)", s.name)

            # Get the last state in each adapter:
            if m is not None and int(m.group(3)) == len(self.adapter_dict[m.group(1)]):
                # Add a transition to the end of the model:
                self.hmm.add_transition(s, self.hmm.end, LibraryModel.SUDDEN_END_PROB)

        # self.dump_as_dotfile()

        # Finalze our HMM:
        self.hmm.bake()

    def dump_as_dotfile(self):
        print("="* 20)
        print("self.hmm")
        print(type(self.hmm))
        print(dir(self.hmm))

        print("="* 20)
        print("self.hmm.graph")
        print(type(self.hmm.graph))
        print(dir(self.hmm.graph))

        print("="* 20)
        print("self.hmm.graph.to_directed()")
        print(type(self.hmm.graph.to_directed()))
        print(dir(self.hmm.graph.to_directed()))

        print("=" * 20)
        # print(self.hmm.to_json())
        # print(networkx.readwrite.gexf.write_gexf(self.hmm.graph, "graph.gexf"))
        # print(networkx.readwrite.graphml.write_graphml(self.hmm.graph, "graph.graphml"))
        # print(networkx.readwrite.gml.write_gml(self.hmm.graph, "graph.gml"))
        # print(networkx.readwrite.json_graph.node_link_data(self.hmm.graph, "graph.json"))

        out_file_name = "graph.dot"
        with open(out_file_name, 'w') as f:
            f.write(f"digraph ")
            f.write(f"{self.name}")
            f.write(" {\n")

            f.write("}")



    def validate_model(self):
        """Ensure that the configuration of this model conforms to the semantics that we require for it to be
        well-formed."""

        # Ensure that our model does not start with a named random segment:
        for name in self.start_element_names:
            if name in self.named_random_segments:
                message = f"ERROR: start segment name is a named random segment: {name}.  " \
                          f"Start segments cannot be named random segments."
                logger.critical(message)
                raise RuntimeError(message)

        # Ensure that our model does not end with a named random segment:
        for name in self.end_element_names:
            if name in self.named_random_segments:
                message = f"ERROR: end segment name is a named random segment: {name}.  " \
                          f"End segments cannot be named random segments."
                logger.critical(message)
                raise RuntimeError(message)

        # Ensure all named random segments have direct connections:
        for name in self.named_random_segments:
            if name not in self.direct_connections_dict:
                message = f"ERROR: Named random segment has no direct connections: {name}.  " \
                          f"All named random segments must be directly connected to another segment."
                logger.critical(message)
                raise RuntimeError(message)

        # Not sure this is needed:
        # # Ensure that all chains of random segments have an "end" at a non-random segment:
        # good_random_nodes = set()
        # for node in self.named_random_segments:
        #     print("="*20)
        #     if self._validate_children_end_in_non_random_nodes(node, good_random_nodes):
        #         good_random_nodes.add(node)
        #     else:
        #         message = f"ERROR: Some named random nodes have connections that do not lead to non-random nodes: " \
        #                   f"{node}"
        #         logger.critical(message)
        #         raise RuntimeError(message)

    # def _validate_children_end_in_non_random_nodes(self, node_name, good_nodes):
    #
    #     print(node_name)
    #
    #     if node_name in good_nodes:
    #         return True
    #     else:
    #         seen_nodes = set()
    #         to_traverse = queue.Queue()
    #         for n in self.direct_connections_dict[node_name]:
    #             to_traverse.put(n)
    #
    #         print(to_traverse.queue)
    #
    #         while not to_traverse.empty():
    #             child = to_traverse.get()
    #             print(f"    {child}")
    #             if child in self.named_random_segments:
    #                 if child in good_nodes:
    #                     print(f"        Good")
    #                 if child in seen_nodes:
    #                     print(f"        Seen")
    #                 else:
    #                     for n in self.direct_connections_dict[child]:
    #                         to_traverse.put(n)
    #             else:
    #                 print(f"        nonrandom")
    #             seen_nodes.add(child)
    #
    #     return False

    @staticmethod
    def _get_state_name(state):
        state_name = state.name

        i = state_name.find(":")
        if i != -1:
            state_name = state_name[:i]

        i = state_name.find("-start")
        if i != -1:
            state_name = state_name[:i]

        i = state_name.find("-end")
        if i != -1:
            state_name = state_name[:i]

        return state_name

    @staticmethod
    def create_model(model_name):
        """Create the model of the given name.
        If the name is not in LibraryModel.has_prebuilt_model, will attempt to read the given model_name as a json
        file to create a model.  This file is assumed to contain the information as per LibraryModel.to_json()."""
        # Get our model:
        if LibraryModel.has_prebuilt_model(model_name):
            logger.info(f"Using %s", LibraryModel.pre_configured_models[model_name]["description"])
            m = LibraryModel.build_pre_configured_model(model_name)
        else:
            logger.info(f"Loading model from json file: %s", model_name)
            m = LibraryModel.from_json_file(model_name)

        return m

    @staticmethod
    def has_prebuilt_model(model_name):
        """Returns True iff LibraryModel.pre_configured_models has a model with the given name.  False otherwise."""
        return model_name in LibraryModel.pre_configured_models

    def _create_key_adapter_order(self):
        """Setup an ordered list of key segments that characterize the correct array element order."""

        # TODO: Generalize this for all library types / segment names!
        # Assumption: self.array_element_structure contains the array elements in our library in the order in which
        #             they appear in the data.
        # Heuristic: The segments that characterize the array elements themselves all have single-character names, so we
        #            filter out all segments from self.array_element_structure with names longer than 1 char.  We then
        #            use these in order to characterize the reads.

        ordered_key_adapters = [s for array in self.array_element_structure
                                for s in array if len(s) == 1]
        return ordered_key_adapters

    @staticmethod
    def _make_fixed_length_random_segment(name, length):
        model = HiddenMarkovModel(name=name)

        state_dict = {}

        random_dist = DiscreteDistribution({
                    "A": LibraryModel.RANDOM_BASE_PROB,
                    "C": LibraryModel.RANDOM_BASE_PROB,
                    "G": LibraryModel.RANDOM_BASE_PROB,
                    "T": LibraryModel.RANDOM_BASE_PROB
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
        model.bake(merge="None")
        return model

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
                        "A": LibraryModel.PER_BASE_MATCH_PROB if target[c] == "A" else LibraryModel.PER_BASE_MISMATCH_PROB,
                        "C": LibraryModel.PER_BASE_MATCH_PROB if target[c] == "C" else LibraryModel.PER_BASE_MISMATCH_PROB,
                        "G": LibraryModel.PER_BASE_MATCH_PROB if target[c] == "G" else LibraryModel.PER_BASE_MISMATCH_PROB,
                        "T": LibraryModel.PER_BASE_MATCH_PROB if target[c] == "T" else LibraryModel.PER_BASE_MISMATCH_PROB,
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
    def _make_random_repeat_model(name):
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

    def to_json(self, outfile=None, indent=4):
        """Serialize this model to a json object and return that json object.
        If outfile is not none, will write the json object to the given file path."""

        model_data = {
            "name": self.name,
            "description": self.description,
            "version": self.version,
            "array_element_structure": self.array_element_structure,
            "adapters": self.adapter_dict,
            "direct_connections": {k: list(v) for k, v in self.direct_connections_dict.items()},
            "start_element_names": list(self.start_element_names),
            "end_element_names": list(self.end_element_names),
            "named_random_segments": list(self.named_random_segments)
        }

        if outfile:
            with open(outfile, 'w') as f:
                json.dump(model_data, f, indent=indent)

        return json.dumps(model_data, indent=indent)

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
            description=json_data["description"],
            version=json_data["version"],
            array_element_structure=tuple(tuple(v) for v in json_data["array_element_structure"]),
            adapters=json_data["adapters"],
            direct_connections={k: set(v) for k, v in json_data["direct_connections"].items()},
            start_element_names=set(json_data["start_element_names"]),
            end_element_names=set(json_data["end_element_names"]),
            named_random_segments=set(json_data["named_random_segments"]),
        )

        return m

    @staticmethod
    def build_pre_configured_model(model_name):
        """Build a pre-configured model based on the given model name.
        If the given name does not appear in self.pre_configured_models, then a KeyError will be thrown."""
        if model_name not in LibraryModel.pre_configured_models:
            raise KeyError(f"Model not found in pre-configured models: {model_name}")

        return LibraryModel(
            name=model_name,
            description=LibraryModel.pre_configured_models[model_name]["description"],
            version=LibraryModel.pre_configured_models[model_name]["version"],
            array_element_structure=LibraryModel.pre_configured_models[model_name]["array_element_structure"],
            adapters=LibraryModel.pre_configured_models[model_name]["adapters"],
            direct_connections=LibraryModel.pre_configured_models[model_name]["direct_connections"],
            start_element_names=LibraryModel.pre_configured_models[model_name]["start_element_names"],
            end_element_names=LibraryModel.pre_configured_models[model_name]["end_element_names"],
            named_random_segments=LibraryModel.pre_configured_models[model_name]["named_random_segments"],
        )

    # TODO: Make an enum for this...
    pre_configured_models = {
        "mas15": {
            "description": "The standard MAS-seq 15 array element model.",
            "version": "2.0.0",
            "array_element_structure": (
                # NOTE: the first element doesn't currently have the "A" adapter in this version of the library.
                ("A", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("B", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("C", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("D", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("E", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("F", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("G", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("H", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("I", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("J", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("K", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("L", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("M", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("N", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                # The last element doesn't currently have the "P" adapter in this version of the library:
                ("O", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter", "P"),
            ),
            "adapters": {
                "10x_Adapter": "TCTACACGACGCTCTTCCGATCT",
                "Poly_A": "A" * 30,
                "3p_Adapter": "GTACTCTGCGTTGATACCACTGCTT",
                "5p_Spacer": "TTTCTTATATGGG",
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
                "CBC": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 16},
                "UMI": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 10},
                "cDNA": RANDOM_SEGMENT_NAME,
            },
            "direct_connections": {
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
                "10x_Adapter": {"CBC"},
                "CBC": {"UMI"},
                "UMI": {"5p_Spacer"},
                "5p_Spacer": {"cDNA"},
                "cDNA": {"Poly_A"},
            },
            "start_element_names": {"A", "10x_Adapter"},
            "end_element_names": {"Poly_A", "P"},
            "named_random_segments": {"UMI", "cDNA", "CBC"},
        },
        "mas15threeP": {
            "description": "The 3' kit MAS-seq 15 array element model.",
            "version": "2.0.0",
            "array_element_structure": (
                # NOTE: the first element doesn't currently have the "A" adapter in this version of the library.
                ("A", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("B", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("C", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("D", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("E", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("F", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("G", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("H", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("I", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("J", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("K", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("L", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("M", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("N", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                # The last element doesn't currently have the "P" adapter in this version of the library:
                ("O", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO", "P"),
            ),
            "adapters": {
                "10x_Adapter": "TCTACACGACGCTCTTCCGATCT",
                "TSO": "CCCATGTACTCTGCGTTGATACCACTGCTT",
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
                "Poly_T": "T" * 30,
                "CBC": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 16},
                "UMI": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 10},
                "cDNA": RANDOM_SEGMENT_NAME,
            },
            "direct_connections": {
                "TSO": {
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
                "10x_Adapter": {"CBC"},
                "CBC": {"UMI"},
                "UMI": {"Poly_T"},
                "Poly_T": {"cDNA"},
                "cDNA": {"TSO"},
            },
            "start_element_names": {"A", "10x_Adapter"},
            "end_element_names": {"P"},
            "named_random_segments": {"UMI", "cDNA", "CBC"},
        },
        "mas10": {
            "description": "The MAS-seq 10 array element model.",
            "version": "2.0.0",
            "array_element_structure": (
                ("Q", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("C", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("M", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("I", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("O", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("J", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("B", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("D", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                ("K", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter"),
                # The last element may not currently have the "R" adapter in this version of the library:
                ("H", "10x_Adapter", "CBC", "UMI", "5p_Spacer", "cDNA", "Poly_A", "3p_Adapter", "R"),
            ),
            "adapters": {
                "10x_Adapter": "TCTACACGACGCTCTTCCGATCT",
                "Poly_A": "A" * 30,
                "3p_Adapter": "GTACTCTGCGTTGATACCACTGCTT",
                "5p_Spacer": "TTTCTTATATGGG",
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
                "CBC": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 16},
                "UMI": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 10},
                "cDNA": RANDOM_SEGMENT_NAME,
            },
            "direct_connections": {
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
                "10x_Adapter": {"CBC"},
                "CBC": {"UMI"},
                "UMI": {"5p_Spacer"},
                "5p_Spacer": {"cDNA"},
                "cDNA": {"Poly_A"},
            },
            "start_element_names": {"Q", "10x_Adapter"},
            "end_element_names": {"Poly_A", "R"},
            "named_random_segments": {"UMI", "cDNA", "CBC"},
        },
        "mas10threeP": {
            "description": "The 3' kit MAS-seq 10 array element model.",
            "version": "2.0.0",
            "array_element_structure": (
                ("Q", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("C", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("M", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("I", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("O", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("J", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("B", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("D", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("K", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                # The last element may not currently have the "R" adapter in this version of the library:
                ("H", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO", "R"),
            ),
            "adapters": {
                "10x_Adapter": "TCTACACGACGCTCTTCCGATCT",
                "TSO": "CCCATGTACTCTGCGTTGATACCACTGCTT",
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
                "Poly_T": "T" * 30,
                "CBC": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 16},
                "UMI": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 10},
                "cDNA": RANDOM_SEGMENT_NAME,
            },
            "direct_connections": {
                "TSO": {
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
                "10x_Adapter": {"CBC"},
                "CBC": {"UMI"},
                "UMI": {"Poly_T"},
                "Poly_T": {"cDNA"},
                "cDNA": {"TSO"},
            },
            "start_element_names": {"Q", "10x_Adapter"},
            "end_element_names": {"R"},
            "named_random_segments": {"UMI", "cDNA", "CBC"},
        },
        "slide-seq": {
            # The slide-seq model is:
            #
            #                 |-----10x_Adapter---->        |--splitter------>               |------Poly_T---------------->                  |--------5p_Adapter----------|                         # noqa
            # AGCTTACTTGTGAAGACTACACGACGCTCTTCCGATCTNNNNNNNNTCTTCAGCGTTCCCGAGANNNNNNNNNNNNNVVTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVNNNNNNNNNNNNNNNNNCCCATGTACTCTGCGTTGATACCACTGCTTACTTGTAAGCTGTCTA...      # noqa
            # |------A------->                      <------|                  <-----------|                                 <----cDNA-------|                              |-------B------>         # noqa
            #                                          V                           V
            #                                    Spatial Barcode 2         Spatial Barcode 1
            "description": "The Slide-seq 15 array element model.",
            "version": "1.0.0",
            "array_element_structure": (
                ("A", "10x_Adapter", "SBC2", "5p_Spacer", "SP1", "UMI", "Poly_T", "cDNA", "5p_Adapter"),
                ("B", "10x_Adapter", "SBC2", "5p_Spacer", "SP1", "UMI", "Poly_T", "cDNA", "5p_Adapter"),
                ("C", "10x_Adapter", "SBC2", "5p_Spacer", "SP1", "UMI", "Poly_T", "cDNA", "5p_Adapter"),
                ("D", "10x_Adapter", "SBC2", "5p_Spacer", "SP1", "UMI", "Poly_T", "cDNA", "5p_Adapter"),
                ("E", "10x_Adapter", "SBC2", "5p_Spacer", "SP1", "UMI", "Poly_T", "cDNA", "5p_Adapter"),
                ("F", "10x_Adapter", "SBC2", "5p_Spacer", "SP1", "UMI", "Poly_T", "cDNA", "5p_Adapter"),
                ("G", "10x_Adapter", "SBC2", "5p_Spacer", "SP1", "UMI", "Poly_T", "cDNA", "5p_Adapter"),
                ("H", "10x_Adapter", "SBC2", "5p_Spacer", "SP1", "UMI", "Poly_T", "cDNA", "5p_Adapter"),
                ("I", "10x_Adapter", "SBC2", "5p_Spacer", "SP1", "UMI", "Poly_T", "cDNA", "5p_Adapter"),
                ("J", "10x_Adapter", "SBC2", "5p_Spacer", "SP1", "UMI", "Poly_T", "cDNA", "5p_Adapter"),
                ("K", "10x_Adapter", "SBC2", "5p_Spacer", "SP1", "UMI", "Poly_T", "cDNA", "5p_Adapter"),
                ("L", "10x_Adapter", "SBC2", "5p_Spacer", "SP1", "UMI", "Poly_T", "cDNA", "5p_Adapter"),
                ("M", "10x_Adapter", "SBC2", "5p_Spacer", "SP1", "UMI", "Poly_T", "cDNA", "5p_Adapter"),
                ("N", "10x_Adapter", "SBC2", "5p_Spacer", "SP1", "UMI", "Poly_T", "cDNA", "5p_Adapter"),
                ("O", "10x_Adapter", "SBC2", "5p_Spacer", "SP1", "UMI", "Poly_T", "cDNA", "5p_Adapter", "P"),
            ),
            "adapters": {
                "10x_Adapter": "TCTACACGACGCTCTTCCGATCT",
                "5p_Adapter": "CCCATGTACTCTGCGTTGATACCACTGCTT",
                "5p_Spacer": "TCTTCAGCGTTCCCGAGA",
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
                "Poly_T": "T" * 30,
                "SBC1": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 6},
                "SBC2": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 8},
                # The UMI might be 7, rather than 9 elements long - not clear from the geneious file.
                "UMI": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 9},
                "cDNA": RANDOM_SEGMENT_NAME,
            },
            "direct_connections": {
                "5p_Adapter": {
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
                "10x_Adapter": {"SBC2"},
                "SBC2": {"5p_Spacer"},
                "5p_Spacer": {"SBC1"},
                "SBC1": {"UMI"},
                "UMI": {"Poly_T"},
                "Poly_T": {"cDNA"},
                "cDNA": {"5p_Adapter"},
            },
            # Right now it won't work properly without both A and 10xAdapter being starts and P and
            # 5pAdapter being ends.
            "start_element_names": {"A", "10x_Adapter"},
            "end_element_names": {"P", "5p_Adapter"},
            "named_random_segments": {"UMI", "SBC2", "SBC1", "cDNA"},
        },
        "mas8prototype": {
            "description": "The prototype MAS-seq 8 array element model.",
            "version": "1.0.0",
            "array_element_structure": (
                # NOTE: the first element may not have the "A" adapter in this version of the library.
                ("A", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("B", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("C", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("D", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("E", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("F", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("G", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO"),
                ("H", "10x_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "TSO", "A"),
            ),
            "adapters": {
                "10x_Adapter": "CTACACGACGCTCTTCCGATCT",
                "Poly_T": "T" * 30,
                "TSO": "CCCATGTACTCTGCGTTGATACCACTGCTT",
                "A": "ACGTACAGT",
                "B": "AAACTGCAT",
                "C": "AGAGTCACT",
                "D": "AGCAAGTTT",
                "E": "AAGTGGTGT",
                "F": "AGTGGACTT",
                "G": "ACAGGTTAT",
                "H": "ATCTCACAT",
                "CBC": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 16},
                "UMI": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 12},
                "cDNA": RANDOM_SEGMENT_NAME,
            },
            "direct_connections": {
                "TSO": {
                    "A",
                    "B",
                    "C",
                    "D",
                    "E",
                    "F",
                    "G",
                    "H",
                },
                "A": {"10x_Adapter"},
                "B": {"10x_Adapter"},
                "C": {"10x_Adapter"},
                "D": {"10x_Adapter"},
                "E": {"10x_Adapter"},
                "F": {"10x_Adapter"},
                "G": {"10x_Adapter"},
                "H": {"10x_Adapter"},
                "10x_Adapter": {"CBC"},
                "CBC": {"UMI"},
                "UMI": {"Poly_T"},
                "Poly_T": {"cDNA"},
                "cDNA": {"TSO"},
            },
            "start_element_names": {"A", "10x_Adapter"},
            "end_element_names": {"TSO", "A"},
            "named_random_segments": {"UMI", "CBC", "cDNA"},
        }
    }


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
