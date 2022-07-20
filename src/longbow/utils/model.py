import re
import logging
import pickle

import click_log

from pomegranate import *

import longbow.utils.constants
from .constants import RANDOM_SEGMENT_NAME, FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME, HPR_SEGMENT_TYPE_NAME, \
    RANDOM_SILENT_STATE_A, RANDOM_SILENT_STATE_B, RANDOM_BASE_STATE, START_STATE_INDICATOR, END_STATE_INDICATOR, \
    BAKE_MERGE_STRATEGY

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger(__name__)
click_log.basic_config(logger)

# # DEBUG:
# logger.setLevel(logging.DEBUG)

starts_with_number_re = re.compile(r"^\d")


class LibraryModel:
    """Model describing a given library preparation.
    The model can annotate the known sections of a read from the library it describes."""

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
                 coding_region,
                 annotation_segments,
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

        self.coding_region = coding_region
        self.annotation_segments = annotation_segments

        self.hmm = None
        self.key_adapters = self._create_key_adapter_order()
        self.key_adapter_set = set(self.key_adapters)

        if do_build:
            self.build()

    @property
    def has_umi_annotation(self):
        return self._has_annotation(longbow.utils.constants.READ_RAW_UMI_TAG) or self._has_annotation(longbow.utils.constants.READ_UMI_TAG)

    @property
    def has_cell_barcode_annotation(self):
        return self._has_annotation(longbow.utils.constants.READ_RAW_BARCODE_TAG) or self._has_annotation(longbow.utils.constants.READ_BARCODE_TAG)

    def _has_annotation(self, annotation_tag):
        if self.annotation_segments:
            # Get all annotation tags:
            for _, tag_tuple_list in self.annotation_segments.items():
                for tag_tuple in tag_tuple_list:
                    if annotation_tag in tag_tuple:
                        return True

        # If we're still running, we didn't find our annotation tag.
        # Therefore we don't have the given annotation.
        return False

    @property
    def has_coding_region(self):
        return self.coding_region is not None

    @property
    def has_named_random_segments(self):
        return self.named_random_segments is not None

    @property
    def num_array_elements(self):
        return len(self.array_element_structure)

    def annotate(self, seq):
        """Annotate the given segment using this model."""
        logp, path = self.hmm.viterbi(seq)

        ppath = []
        for p, (idx, state) in enumerate(path[1:-1]):
            if (
                    not state.name.endswith(START_STATE_INDICATOR)
                    and not state.name.endswith(END_STATE_INDICATOR)
                    and ":RD" not in state.name
                    and ":D" not in state.name
            ):
                ppath.append(f"{re.split(':', state.name)[0]}")

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
                    # TODO: Replace all references to "5p_Adapter" with references to the models themselves.
                    if allow_missing_first_adapter and (key_adapter_indx == 1) and \
                            (ordered_segment_names[0] == "5p_Adapter"):
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

    def get_all_node_names(self):
        """Return a list containing the names of all nodes in this model."""

        # Let's be exhaustive here just to make sure we don't miss anything:
        all_node_names = {n for element in self.array_element_structure for n in element}

        if self.named_random_segments is not None:
            for n in self.named_random_segments:
                all_node_names.add(n)

        for n in self.start_element_names:
            all_node_names.add(n)

        for n in self.end_element_names:
            all_node_names.add(n)

        for k in self.direct_connections_dict:
            all_node_names.add(k)
            for n in self.direct_connections_dict[k]:
                all_node_names.add(n)

        return sorted(list(all_node_names))

    def build(self):
        """Build the HMM underlying this model given our segment information."""

        # Warn if model is deprecated
        try:
            if self.pre_configured_models[self.name]['deprecated']:
                logger.warning(f"Model {self.name} is deprecated.")
        except KeyError:
            # This is OK - it may be the user specifying a model that is not in our preconfigured set.
            pass

        # Validate our model here so we can go through and just worry about creating it:
        self.validate_model()

        # Initialize the states in our model:
        self._initialize_hmm_states()

        # Initialize the transitions in our model:
        # self._initialize_hmm_transitions()
        self._initialize_hmm_transitions_new()

        # DEBUGGING:
        if logger.isEnabledFor(logging.DEBUG):
            self.dump_as_dotfile(do_subgraphs=False)
            self.dump_as_dotfile_simple()
            self.to_json(f"longbow_model_{self.name}.v{self.version}.json")

            with open(f"longbow_model_{self.name}.v{self.version}.dense_transition_matrix.pickle", 'wb') as f:
                pickle.dump(self.hmm.dense_transition_matrix(), f)
            with open(f"longbow_model_{self.name}.v{self.version}.emission_distributions.txt", 'w') as f:
                print(self.hmm, file=f, flush=True)

        # Finalze our HMM:
        self.hmm.bake(merge=BAKE_MERGE_STRATEGY)

    def _get_segment_start_end_rda_and_rdb_nodes(self):
        """Create a Dictionary of model segment starting states, random start, and random end"""
        segment_starts = {}
        segment_ends = {}
        rda = None
        rdb = None
        for s in self.hmm.states:
            if self._is_state_last_in_section(s):
                segment_ends[s.name] = s
            if s.name.endswith(START_STATE_INDICATOR) and not s.name.startswith(RANDOM_SEGMENT_NAME):
                segment_starts[re.sub(START_STATE_INDICATOR, "", s.name)] = s
            elif s.name.startswith(f"{RANDOM_SEGMENT_NAME}:{RANDOM_SILENT_STATE_A}"):
                rda = s
            elif s.name.startswith(f"{RANDOM_SEGMENT_NAME}:{RANDOM_SILENT_STATE_B}"):
                rdb = s

        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"Segment starts: {[s.name for s in segment_starts.values()]}")
            logger.debug(f"RDA: {rda}")
            logger.debug(f"RDB: {rdb}")

        return segment_starts, segment_ends, rda, rdb

    def _create_random_to_segment_start_transitions(self, segment_starts, rda):
        """Create transitions from rda to all segments' start nodes.
        NOTE: transitions are NOT made from rda to named random segments."""

        # link random element model silent state A to all of our starts
        # This is the transition from `random` to one of our array elements.
        # This only happens for models with named random segments:
        per_state_transition_prob = 1.0
        if self.has_named_random_segments:
            # Get the number of "standard" segments.
            # i.e. states that are not named random segments:
            num_standard_states = 0
            for sname in segment_starts:
                if sname not in self.named_random_segments:
                    num_standard_states += 1
            per_state_transition_prob = 1.0 / num_standard_states

        for sname in segment_starts:
            # We do NOT account for rda -> start_element transitions here.
            # Those should be accounted for in another method.
            #
            # We also exclude named random segments because they model the same thing.
            if sname not in self.start_element_names and \
                    (self.has_named_random_segments and sname not in self.named_random_segments):
                if logger.isEnabledFor(logging.DEBUG):
                    logger.debug(
                        f"Adding transition: {self.hmm.start.name} -> {segment_starts[sname].name} "
                        f"@ {per_state_transition_prob:.3f}")
                self.hmm.add_transition(rda, segment_starts[sname], per_state_transition_prob)

    def _create_array_start_and_end_transitions(self, segment_starts, rda):
        """Create transitions between:
        self.model.start -> [ALL starting segments' start nodes]
        rda -> [ALL ending segments' start nodes]"""

        # Link the model's silent start state to all of our library start nodes.
        # (These are not the same as start nodes for each segment.  This is the start of the array.)
        # NOTE: we cannot start in a named random segment.
        for sname in segment_starts:
            if sname in self.start_element_names:
                weight = LibraryModel.START_RANDOM_PROB / len(self.start_element_names)
                if logger.isEnabledFor(logging.DEBUG):
                    logger.debug(
                        f"Adding transition: {self.hmm.start.name} -> {segment_starts[sname].name} @ {weight:.3f}")
                self.hmm.add_transition(self.hmm.start, segment_starts[sname], weight)

        # Link random element model silent state A to all of our end elements.
        # This is the transition from `random` to an element that ends an array.
        # (i.e. these are not the same as end nodes for each segment.  This is the end of the array.)
        # NOTE: we cannot end in a named random segment.
        for sname in segment_starts:
            if sname in self.end_element_names:
                if logger.isEnabledFor(logging.DEBUG):
                    logger.debug(
                        f"Adding transition: RDA -> {segment_starts[sname].name} @ {1.0 / len(self.end_element_names):.3f}")
                self.hmm.add_transition(rda, segment_starts[sname], 1.0 / len(self.end_element_names))

    def _initialize_hmm_transitions(self):
        logger.debug("Identifying actionable nodes from model...")
        segment_starts, segment_ends, rda, rdb = self._get_segment_start_end_rda_and_rdb_nodes()

        # Create transitions for the array start and array end:
        logger.debug("Creating start/end library model transitions...")
        self._create_array_start_and_end_transitions(segment_starts, rda)

        # link up adapter final states according to our direct connections dictionary
        logger.debug("Adding transitions for direct connections...")
        for _, s in segment_ends.items():
            logger.debug("State: %s", s.name)
            sname = LibraryModel._get_state_base_name(s)

            # If we have direct connections from this state to some other states then we add
            # them into the model:
            if sname in self.direct_connections_dict:
                for dcname in self.direct_connections_dict[sname]:
                    if logger.isEnabledFor(logging.DEBUG):
                        logger.debug(
                            f"    Adding transition: {s.name} -> {segment_starts[dcname].name} @ {(1.0 - LibraryModel.SUDDEN_END_PROB) / len(self.direct_connections_dict[sname]):.3f}")
                    self.hmm.add_transition(
                        s, segment_starts[dcname],
                        1.0 / len(self.direct_connections_dict[sname])
                    )

            # If we have no direct connections from this state, we add in a transition from this
            # state to the silent random B state:
            else:
                # Verify this probability is the notional equivalent:
                if logger.isEnabledFor(logging.DEBUG):
                    logger.debug(
                        f"    Adding transition (no DCs): {s.name} -> RDB @ {LibraryModel.RAND_RAND_PROB:.3f}")
                self.hmm.add_transition(s, rdb, LibraryModel.RAND_RAND_PROB)

        logger.debug("Adding transitions for short-circuit terminations...")
        for _, s in segment_ends.items():
            logger.debug("State: %s", s.name)
            self._create_sudden_end_transitions(s)

    def _initialize_hmm_transitions_new(self):
        logger.debug("Identifying actionable nodes from model...")
        segment_starts, segment_ends, rda, rdb = self._get_segment_start_end_rda_and_rdb_nodes()

        # Create transitions for the array start and array end:
        logger.debug("Creating start/end library model transitions...")
        self._create_array_start_and_end_transitions(segment_starts, rda)

        # Create transitions from RDA to segment starts:
        # NOTE: This will NOT create transitions to start nodes of named random segments.
        logger.debug("Creating random -> segment start node transitions...")
        self._create_random_to_segment_start_transitions(segment_starts, rda)

        # link up adapter final states according to our direct connections dictionary
        logger.debug("Adding transitions for direct connections...")
        for _, s in segment_ends.items():
            logger.debug("State: %s", s.name)
            sname = LibraryModel._get_state_base_name(s)

            # If we have direct connections from this state to some other states then we add
            # them into the model:
            if sname in self.direct_connections_dict:
                for dcname in self.direct_connections_dict[sname]:
                    if logger.isEnabledFor(logging.DEBUG):
                        logger.debug(
                            f"    Adding transition: {s.name} -> {segment_starts[dcname].name} @ {(1.0 - LibraryModel.SUDDEN_END_PROB) / len(self.direct_connections_dict[sname]):.3f}")
                    self.hmm.add_transition(
                        s, segment_starts[dcname],
                        # Original transition prob:
                        1.0 / len(self.direct_connections_dict[sname])
                        # Attempted new transition prob:
                        # (1.0 - LibraryModel.SUDDEN_END_PROB) / len(self.direct_connections_dict[sname])
                    )

                # # Also add in a connection to the silent random B state if this is a non-random section
                # # and it does not connect to a named random section:
                # has_direct_connection_to_named_random = False
                # if self.has_named_random_segments and (sname not in self.named_random_segments):
                #     for n in self.direct_connections_dict[sname]:
                #         if n in self.named_random_segments:
                #             has_direct_connection_to_named_random = True
                #             break
                # if not has_direct_connection_to_named_random:
                #     # Verify this probability is the notional equivalent:
                #     if logger.isEnabledFor(logging.DEBUG):
                #         logger.debug(
                #             f"    Adding transition (no rand connections): {s.name} -> RDB @ {LibraryModel.RAND_RAND_PROB:.3f}")
                #     self.hmm.add_transition(s, rdb, LibraryModel.RAND_RAND_PROB)

            # If we have no direct connections from this state, we add in a transition from this
            # state to the silent random B state:
            else:
                # Verify this probability is the notional equivalent:
                if logger.isEnabledFor(logging.DEBUG):
                    logger.debug(
                        f"    Adding transition (no DCs): {s.name} -> RDB @ {LibraryModel.RAND_RAND_PROB:.3f}")
                self.hmm.add_transition(s, rdb, LibraryModel.RAND_RAND_PROB)

        logger.debug("Adding transitions for short-circuit terminations...")
        for _, s in segment_ends.items():
            logger.debug("State: %s", s.name)
            self._create_sudden_end_transitions(s)

    def _create_sudden_end_transitions(self, s):
        # Link up all adapters to model end state;
        # Add a transition to the end of the model:
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(
                f"    Adding transition: {s.name} -> {self.hmm.end.name} @ {LibraryModel.SUDDEN_END_PROB:.3f}")
        self.hmm.add_transition(s, self.hmm.end, LibraryModel.SUDDEN_END_PROB)

    def _initialize_hmm_states(self):
        # Create the centralized random repeat model that allows for errors between the known markers:
        self.hmm = LibraryModel._make_random_repeat_model()
        # Create models for individual adapter sequences,
        # taking special care with named random sections:
        for k, v in self.adapter_dict.items():
            adapter_hmm = None
            if type(v) is str:
                if self.has_named_random_segments and (k in self.named_random_segments):
                    # We have a named random segment.
                    # We need to make a named random model:
                    adapter_hmm = LibraryModel._make_named_random_model(name=k)
                else:
                    # This must be a normal string of bases to model:
                    adapter_hmm = LibraryModel._make_global_alignment_model(v, k)
            elif type(v) is dict:
                segment_type = list(v.keys())[0]

                # NOTE: Must add in new model types here when we add them:
                if segment_type == FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME:
                    # Fixed length random segments for barcodes (and similar):
                    adapter_hmm = LibraryModel._make_fixed_length_random_segment(k, list(v.values())[0])
                elif segment_type == HPR_SEGMENT_TYPE_NAME:
                    # Homopolymer repeat region:
                    base, hpr_length = list(v.values())[0]
                    adapter_hmm = LibraryModel._make_homopolymer_repeat_model(k, base, hpr_length)
                else:
                    logger.critical(f"Unknown special model type: {segment_type}")
                    raise RuntimeError(f"Unknown special model type: {segment_type}")
            self.hmm.add_model(adapter_hmm)
        # Finalize all models that we have for now:
        self.hmm.bake(merge=BAKE_MERGE_STRATEGY)

    def _is_state_last_in_section(self, s):
        """Returns True if the given state is considered to be the last state in its segment."""

        # This is the old logic for determining if a state is the last one in a segment:
        if s.name == self.hmm.end.name:
            return False

        # Get the base name of the node:
        base_name = LibraryModel._get_state_base_name(s.name)

        # Decode the segment number from the name:
        m = re.match(r"^(\w+):([MID])(\d+)", s.name)
        if m is None and (base_name == RANDOM_SEGMENT_NAME or
                          (self.has_named_random_segments and base_name not in self.named_random_segments)):
            return False
        # Allow for the state number to be 0 in the case of a named random segment (e.g. cDNA):
        segment_state_num = int(m.group(3)) if m is not None else 0

        # Case for normal segments:
        if type(self.adapter_dict[base_name]) == str:
            if self.adapter_dict[base_name] == RANDOM_SEGMENT_NAME \
                    and self.has_named_random_segments \
                    and base_name in self.named_random_segments:

                state_end_name = f"{base_name}{END_STATE_INDICATOR}"

                # The final state for the named random segment types is the end state:
                return s.name == state_end_name
            else:
                # Check to see if it's the last state in the segment:
                return (m is not None) and (segment_state_num == len(self.adapter_dict[base_name]))
        # case for special segment types (all of which are dictionaries):
        else:
            segment_type = list(self.adapter_dict[base_name].keys())[0]
            if segment_type == FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME:
                # For fixed length random segments, the value for the type dict is the length:
                segment_length = self.adapter_dict[base_name][segment_type]
            elif segment_type == HPR_SEGMENT_TYPE_NAME:
                # For homopolymer repeat segments, the value for the type dict is a tuple of (NUCLEOTIDE, LENGTH):
                segment_length = self.adapter_dict[base_name][segment_type][1]
            else:
                logger.error(f"UNKNOWN TYPE OF SEGMENT: {segment_type}")
                return False

            return segment_state_num == segment_length

        # # New logic for determining the above:
        # # The only restriction on reporting an end state is that we don't allow the final state of the model to
        # # be an end state.
        # return s.name.endswith(END_STATE_INDICATOR) and (s.name != self.hmm.end.name)

    def dump_as_dotfile_simple(self, out_file_name=None):
        """Dump this LibraryModel to disk as a simple dot file for visualization."""

        if out_file_name is None:
            out_file_name = f"longbow_model_{self.name}.v{self.version}.simple.dot"

        # Get the model as a dictionary so we can handle it
        model_dict = self.hmm.to_dict()

        # Save a map here of node names so we can make some subgraphs:
        named_node_state_maps = {n: [] for n in self.get_all_node_names()}
        named_node_state_maps[RANDOM_SEGMENT_NAME] = []

        with open(out_file_name, 'w') as f:
            f.write(f"digraph ")
            f.write(f"{self.name}")
            f.write(" {\n")

            f.write("\n")

            # Write out the nodes:
            f.write("    // States:\n")
            states = set()
            for node in model_dict['states']:
                name = LibraryModel._get_state_base_name(node['name'])
                states.add(name)
            for state in states:
                f.write(f"    \"{state}\";\n")
            f.write("\n")

            # Write out the edges:
            f.write("    // Transitions:\n")
            transitions = dict()
            for edge in model_dict['edges']:
                source = model_dict['states'][edge[0]]
                dest = model_dict['states'][edge[1]]
                weight = edge[2]

                source_name = LibraryModel._get_state_base_name(source['name'])
                dest_name = LibraryModel._get_state_base_name(dest['name'])

                # TODO: This may be too aggressive a filtering criterion.
                if source_name == dest_name:
                    continue

                transition_descriptor = f"\"{source_name}\" -> \"{dest_name}\""

                transitions[transition_descriptor] =weight
            for t, w in transitions.items():
                f.write(f"    {t} [weight={w:.3f}];\n")

            f.write("}\n")

    def dump_as_dotfile(self, out_file_name=None, do_subgraphs=True):
        """Dump this LibraryModel to disk as a dot file for visualization."""

        if out_file_name is None:
            out_file_name = f"longbow_model_{self.name}.v{self.version}.dot"

        # Get the model as a dictionary so we can handle it
        model_dict = self.hmm.to_dict()

        # Save a map here of node names so we can make some subgraphs:
        named_node_state_maps = {n: [] for n in self.get_all_node_names()}
        named_node_state_maps[RANDOM_SEGMENT_NAME] = []

        with open(out_file_name, 'w') as f:
            f.write(f"digraph ")
            f.write(f"{self.name}")
            f.write(" {\n")

            f.write("\n")

            # Write out the nodes:
            f.write("    // States:\n")
            if do_subgraphs:
                # Collect our subgraphs:
                for node in model_dict['states']:
                    name = LibraryModel._get_state_base_name(node['name'])
                    named_node_state_maps[name].append(node['name'])

                # Write out the subgraphs we have found:
                for subgraph_name, node_names in named_node_state_maps.items():
                    f.write(f"    subgraph \"{subgraph_name}\"")
                    f.write(" {\n")
                    for node_name in node_names:
                        f.write(f"        \"{node_name}\";\n")
                    f.write("    }\n")
            else:
                for node in model_dict['states']:
                    f.write(f"    \"{node['name']}\";\n")

            f.write("\n")

            # Write out the edges:
            f.write("    // Transitions:\n")
            for edge in model_dict['edges']:
                source = model_dict['states'][edge[0]]
                dest = model_dict['states'][edge[1]]
                weight = edge[2]

                f.write(f"    \"{source['name']}\" -> \"{dest['name']}\" [weight={weight:.3f}];\n")

            f.write("}\n")

    def validate_model(self):
        """Ensure that the configuration of this model is well-formed and conforms to the semantics that we require."""

        # Ensure that our model does not start with a named random segment:
        for name in self.start_element_names:
            if self.has_named_random_segments and (name in self.named_random_segments):
                message = f"ERROR: start segment name is a named random segment: {name}.  " \
                          f"Start segments cannot be named random segments."
                logger.critical(message)
                raise RuntimeError(message)

        # Ensure that our model does not end with a named random segment:
        for name in self.end_element_names:
            if self.has_named_random_segments and (name in self.named_random_segments):
                message = f"ERROR: end segment name is a named random segment: {name}.  " \
                          f"End segments cannot be named random segments."
                logger.critical(message)
                raise RuntimeError(message)

        # Ensure all named random segments have direct connections:
        if self.has_named_random_segments:
            for name in self.named_random_segments:
                if name not in self.direct_connections_dict:
                    message = f"ERROR: Named random segment has no direct connections: {name}.  " \
                              f"All named random segments must be directly connected to another segment."
                    logger.critical(message)
                    raise RuntimeError(message)

        # Ensure that the coding region is specified somewhere in the model as a string:
        if self.coding_region is not None:
            if type(self.coding_region) is not str:
                message = f"ERROR: coding_region must be a string."
                logger.critical(message)
                raise RuntimeError(message)
            elif self.coding_region not in self.adapter_dict:
                message = f"ERROR: coding_region name does not appear in adapter dictionary: " \
                          f"{self.coding_region} not in ({self.adapter_dict.keys()})."
                logger.critical(message)
                raise RuntimeError(message)

        # Ensure that all segments in the array element structure show up in the adapter dictionary:
        for array_element in self.array_element_structure:
            for adapter in array_element:
                # `random` nodes are OK in the array element structure:
                if adapter not in self.adapter_dict and adapter != RANDOM_SEGMENT_NAME:
                    message = f"ERROR: adapter from array element structure does not appear in adapter dictionary: " \
                              f"{adapter} not in ({self.adapter_dict.keys()})."
                    logger.critical(message)
                    raise RuntimeError(message)

        # Ensure that all segments in the direct connections dict exist in the adapter dictionary:
        for source, destination_nodes in self.direct_connections_dict.items():
            if source not in self.adapter_dict:
                message = f"ERROR: source adapter from direct connections does not appear in adapter dictionary: " \
                          f"{source} not in ({self.adapter_dict.keys()})."
                logger.critical(message)
                raise RuntimeError(message)
            for dest in destination_nodes:
                if dest not in self.adapter_dict:
                    message = f"ERROR: dest adapter from direct connections does not appear in adapter dictionary: " \
                              f"{dest} not in ({self.adapter_dict.keys()})."
                    logger.critical(message)
                    raise RuntimeError(message)

    @staticmethod
    def _get_state_base_name(state):
        if type(state) is str:
            state_name = state
        else:
            state_name = state.name

        i = state_name.find(":")
        if i != -1:
            state_name = state_name[:i]

        i = state_name.find(START_STATE_INDICATOR)
        if i != -1:
            state_name = state_name[:i]

        i = state_name.find(END_STATE_INDICATOR)
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
    def _make_homopolymer_repeat_model(name, nucleotide, expected_length):
        # TODO: Do a better bake-off for the HPR models.
        #  Currently v1 works very well, but you should still verify it's the best.
        #  Actually, it should be refined to "center" on the specified length and have "extra"
        #  transitions at the 90% and 110% marks.  This is similar to what you're doing now,
        #  but a little different.

        # return LibraryModel._make_homopolymer_repeat_model_v0(name, nucleotide, expected_length)
        return LibraryModel._make_homopolymer_repeat_model_v1(name, nucleotide, expected_length)
        # return LibraryModel._make_homopolymer_repeat_model_v2(name, nucleotide, expected_length)

    @staticmethod
    def _make_homopolymer_repeat_model_v0(name, nucleotide, expected_length):
        logger.debug("Making Model: HOMOPOLYMER_REPEAT (%s:%s x %d)", name, nucleotide, expected_length)

        # The homopolymer repeat model begins its life as a global alignment model consisting
        # of one base repeated for the expected length of the HPR:
        model = LibraryModel._make_global_alignment_model(nucleotide * expected_length, name)

        # Now we add extra transitions to help with HPR detection:
        for s in model.states:
            # Add a transition from each match state to the silent end state to allow for a
            # faster short-circuit end to the HPR that's less than the expected length:
            if s.name.startswith(f"{name}:M"):
                # We don't need an extra transition to the end here:
                if s.name != f"{name}:M{expected_length}":
                    model.add_transition(s, model.end, LibraryModel.HPR_SUDDEN_END_PROB)

                # Add a transition from each match state back to the first Match state to more easily allow for
                # longer than expected HPRs:
                model.add_transition(s, model.start, LibraryModel.HPR_MODEL_RECURRENCE_PROB)

            # Add a transition from each insert state back to the first Match state to more easily allow for
            # longer than expected HPRs:
            elif s.name.startswith(f"{name}:I"):
                # We don't need a recurrence transition at insert state 0
                # (arbitrary, but it makes sense intuitively to me).
                if s.name != f"{name}:I0":
                    model.add_transition(s, model.start, LibraryModel.HPR_MODEL_RECURRENCE_PROB)
                if s.name != f"{name}:I{expected_length}":
                    model.add_transition(s, model.end, LibraryModel.HPR_SUDDEN_END_PROB)

        # Finalize and return the model:
        model.bake(merge=BAKE_MERGE_STRATEGY)
        return model

    @staticmethod
    def _make_homopolymer_repeat_model_v1(name, nucleotide, expected_length):
        logger.debug("Making Model: HOMOPOLYMER_REPEAT (%s:%s x %d)", name, nucleotide, expected_length)

        model = HiddenMarkovModel(name=name)
        s = {}

        # Define our distributions here so we can use them later:
        random_distribution = DiscreteDistribution({
                "A": LibraryModel.RANDOM_BASE_PROB,
                "C": LibraryModel.RANDOM_BASE_PROB,
                "G": LibraryModel.RANDOM_BASE_PROB,
                "T": LibraryModel.RANDOM_BASE_PROB
            })

        standard_distribution = DiscreteDistribution(
                {
                    "A": LibraryModel.PER_BASE_MATCH_PROB if nucleotide == "A" else LibraryModel.PER_BASE_MISMATCH_PROB,
                    "C": LibraryModel.PER_BASE_MATCH_PROB if nucleotide == "C" else LibraryModel.PER_BASE_MISMATCH_PROB,
                    "G": LibraryModel.PER_BASE_MATCH_PROB if nucleotide == "G" else LibraryModel.PER_BASE_MISMATCH_PROB,
                    "T": LibraryModel.PER_BASE_MATCH_PROB if nucleotide == "T" else LibraryModel.PER_BASE_MISMATCH_PROB,
                }
            )

        # Add an extra penalty to the last state for being a base that isn't the one in this HPR:
        bookend_state_distribution = DiscreteDistribution(
                {
                    "A": LibraryModel.HPR_BOOKEND_MATCH_PROB if nucleotide == "A" else LibraryModel.HPR_BOOKEND_MISMATCH_PROB,
                    "C": LibraryModel.HPR_BOOKEND_MATCH_PROB if nucleotide == "C" else LibraryModel.HPR_BOOKEND_MISMATCH_PROB,
                    "G": LibraryModel.HPR_BOOKEND_MATCH_PROB if nucleotide == "G" else LibraryModel.HPR_BOOKEND_MISMATCH_PROB,
                    "T": LibraryModel.HPR_BOOKEND_MATCH_PROB if nucleotide == "T" else LibraryModel.HPR_BOOKEND_MISMATCH_PROB,
                }
            )

        # add states
        i0 = State(random_distribution, name=f"{name}:I0")

        model.add_state(i0)
        s[i0.name] = i0

        for c in range(expected_length):
            if c ==0 or c == expected_length-1:
                # mc = State(standard_distribution, name=f"{name}:M{c + 1}")
                mc = State(bookend_state_distribution, name=f"{name}:M{c + 1}")
            else:
                mc = State(standard_distribution, name=f"{name}:M{c + 1}")

            ic = State(random_distribution, name=f"{name}:I{c + 1}")

            model.add_states([mc, ic])

            s[mc.name] = mc
            s[ic.name] = ic

        # Add transitions for starting states:
        model.add_transition(model.start, s[f"{name}:I0"], LibraryModel.MATCH_INDEL_PROB)
        model.add_transition(model.start, s[f"{name}:M1"], LibraryModel.MATCH_MATCH_PROB)
        model.add_transition(s[f"{name}:M1"], model.end, LibraryModel.HPR_SUDDEN_END_PROB)
        model.add_transition(s[f"{name}:M1"], model.start, LibraryModel.HPR_MODEL_RECURRENCE_PROB)

        model.add_transition(s[f"{name}:I0"], s[f"{name}:I0"], LibraryModel.INDEL_CONTINUATION_PROB)
        model.add_transition(s[f"{name}:I0"], s[f"{name}:M1"], LibraryModel.INDEL_SWITCH_PROB / 2)

        # Add transitions for middle states:
        for c in range(1, expected_length):
            model.add_transition(s[f"{name}:I{c}"], s[f"{name}:I{c}"], LibraryModel.INDEL_SWITCH_PROB / 2)
            model.add_transition(s[f"{name}:I{c}"], s[f"{name}:M{c + 1}"], LibraryModel.INDEL_CONTINUATION_PROB)
            model.add_transition(s[f"{name}:I{c}"], model.start, LibraryModel.HPR_MODEL_RECURRENCE_PROB)
            model.add_transition(s[f"{name}:I{c}"], model.end, LibraryModel.HPR_SUDDEN_END_PROB)

            model.add_transition(s[f"{name}:M{c}"], s[f"{name}:I{c}"], LibraryModel.MATCH_INDEL_PROB)
            model.add_transition(s[f"{name}:M{c}"], s[f"{name}:M{c + 1}"], LibraryModel.MATCH_MATCH_PROB)

        # Add transitions for the last states:
        model.add_transition(s[f"{name}:I{expected_length}"], s[f"{name}:I{expected_length}"],
                             LibraryModel.INDEL_SWITCH_PROB / 2)
        model.add_transition(
            s[f"{name}:I{expected_length}"],
            model.end,
            LibraryModel.INDEL_CONTINUATION_PROB + (LibraryModel.INDEL_SWITCH_PROB / 2)
        )

        model.add_transition(
            s[f"{name}:M{expected_length}"], s[f"{name}:I{expected_length}"], LibraryModel.MATCH_TRAIL_INSERT_PROB
        )
        model.add_transition(s[f"{name}:M{expected_length}"], model.end, LibraryModel.MATCH_END_PROB)
        model.add_transition(s[f"{name}:M{expected_length}"], model.start, LibraryModel.MATCH_END_PROB)

        # Finalize and return the model:
        model.bake(merge=BAKE_MERGE_STRATEGY)
        return model

    @staticmethod
    def _make_homopolymer_repeat_model_v2(name, nucleotide, expected_length):
        logger.debug("Making Model: HOMOPOLYMER_REPEAT (%s:%s x %d)", name, nucleotide, expected_length)

        model = HiddenMarkovModel(name=name)
        states = dict()

        # Add states and "internal" transitions:
        for i in range(1, expected_length+1):

            state_name = f"{name}:I{i}"

            s = State(
                DiscreteDistribution(
                    {
                        "A": LibraryModel.HPR_MATCH_PROB if nucleotide == "A" else LibraryModel.HPR_MISMATCH_PROB,
                        "C": LibraryModel.HPR_MATCH_PROB if nucleotide == "C" else LibraryModel.HPR_MISMATCH_PROB,
                        "G": LibraryModel.HPR_MATCH_PROB if nucleotide == "G" else LibraryModel.HPR_MISMATCH_PROB,
                        "T": LibraryModel.HPR_MATCH_PROB if nucleotide == "T" else LibraryModel.HPR_MISMATCH_PROB,
                    }
                ),
                name=f"{name}:I{i}",
            )

            model.add_state(s)
            states[state_name] = s

            # Add transitions with probabilities based on whether this is a "short-circuit" end node or not:
            if i == expected_length:
                end_transition_prob = 0.1
            else:
                end_transition_prob = 0.01

            model.add_transition(s, model.end, end_transition_prob)

            # Add self-transition:
            model.add_transition(s, s, (1-end_transition_prob)/2)

            # Add transitions between consecutive states:
            if i > 1:
                model.add_transition(states[f"{name}:I{i - 1}"], s, (1-end_transition_prob)/2)

        # Add start transition:
        model.add_transition(model.start, states[f"{name}:I1"], 1)

        # Finalize and return the model:
        model.bake(merge=BAKE_MERGE_STRATEGY)
        return model

    @staticmethod
    def _make_fixed_length_random_segment(name, length):
        logger.debug("Making Model: FIXED_LENGTH_RANDOM (%s:%d)", name, length)
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
        model.bake(merge=BAKE_MERGE_STRATEGY)
        return model

    @staticmethod
    def _make_global_alignment_model(target, name=None):
        logger.debug("Making Model: GLOBAL_ALIGNMENT (%s)", name)
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

        model.bake(merge=BAKE_MERGE_STRATEGY)
        return model

    @staticmethod
    def _make_named_random_model(name):
        logger.debug("Making Model: NAMED RANDOM (%s)", name)

        model = HiddenMarkovModel(name=name)

        # This is a VERY simple repeat model.
        # The idea here is that is a user specifies a named random segment, then we know
        # there should be at least one random base here.
        ri = State(
            DiscreteDistribution({
                "A": LibraryModel.RANDOM_BASE_PROB,
                "C": LibraryModel.RANDOM_BASE_PROB,
                "G": LibraryModel.RANDOM_BASE_PROB,
                "T": LibraryModel.RANDOM_BASE_PROB
            }),
            name=f"{name}:{RANDOM_BASE_STATE}",
        )

        model.add_state(ri)

        model.add_transition(model.start, ri, 1)
        model.add_transition(ri, ri, LibraryModel.NAMED_RAND_CONTINUE_PROB)
        model.add_transition(ri, model.end, LibraryModel.NAMED_RAND_EXIT_PROB)

        model.bake(merge=BAKE_MERGE_STRATEGY)
        return model

    @staticmethod
    def _make_random_repeat_model():
        logger.debug(f"Making Model: RANDOM REPEAT")
        model = HiddenMarkovModel(name=RANDOM_SEGMENT_NAME)

        # add states
        ri = State(
            DiscreteDistribution({
                "A": LibraryModel.RANDOM_BASE_PROB,
                "C": LibraryModel.RANDOM_BASE_PROB,
                "G": LibraryModel.RANDOM_BASE_PROB,
                "T": LibraryModel.RANDOM_BASE_PROB
            }),
            name=f"{RANDOM_SEGMENT_NAME}:{RANDOM_BASE_STATE}",
        )
        rda = State(None, name=f"{RANDOM_SEGMENT_NAME}:{RANDOM_SILENT_STATE_A}")
        rdb = State(None, name=f"{RANDOM_SEGMENT_NAME}:{RANDOM_SILENT_STATE_B}")

        model.add_states([ri, rda, rdb])

        # add transitions
        model.add_transition(model.start, rda, LibraryModel.START_AND_END_RANDOM_PROB)
        model.add_transition(model.start, ri, LibraryModel.START_AND_END_RANDOM_PROB)

        model.add_transition(ri, ri, LibraryModel.RAND_INS_CONTINUATION_PROB)
        model.add_transition(ri, rda, LibraryModel.RAND_INS_TO_DEL_PROB)
        model.add_transition(ri, model.end, LibraryModel.RAND_INS_END_PROB)

        model.add_transition(rdb, ri, LibraryModel.RAND_RAND_PROB)
        model.add_transition(rdb, model.end, LibraryModel.START_AND_END_RANDOM_PROB)

        model.bake(merge=BAKE_MERGE_STRATEGY)

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
        }

        # Set our new "optional" fields here:
        if self.has_named_random_segments:
            model_data["named_random_segments"] = list(self.named_random_segments)

        if self.coding_region is not None:
            model_data["coding_region"] = self.coding_region

        if self.annotation_segments is not None:
            model_data["annotation_segments"] = {k: list(v) for k, v in self.annotation_segments.items()}

        if outfile:
            with open(outfile, 'w') as f:
                json.dump(model_data, f, indent=indent)

        return json.dumps(model_data, indent=indent)

    def has_annotation_tag(self, tag):
        """Returns true if this model has an annotation segment that is written out to the given tag."""
        return self.get_segment_name_for_annotation_tag(tag) is not None

    def get_segment_name_for_annotation_tag(self, tag):
        """Returns the name of the model segment that corresponds to the given tag if such a segment exists.
        Otherwise returns None."""
        barcode_seg_name = None
        for n, tag_tuple_list in self.annotation_segments.items():
            for seg_tag, seg_pos_tag in tag_tuple_list:
                if seg_tag == tag:
                    barcode_seg_name = n

            if barcode_seg_name:
                break

        return barcode_seg_name

    def get_segment_length(self, segment_name):
        """Returns the number of bases in this model for the given segment_name if that segment_name is in this model.
        If the segment_name is a random segment with variable length, zero is returned.
        If the segment_name is not in this model, None is returned."""
        if segment_name not in self.adapter_dict.keys():
            return None

        segment_length = None
        if type(self.adapter_dict[segment_name]) == str:
            segment_length = len(self.adapter_dict[segment_name])
        elif type(self.adapter_dict[segment_name]) == dict:
            seg_type = next(iter(self.adapter_dict[segment_name].keys()))
            if seg_type == longbow.utils.constants.HPR_SEGMENT_TYPE_NAME:
                segment_length = next(iter(self.adapter_dict[segment_name].values()))[1]
            elif seg_type == longbow.utils.constants.FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME:
                segment_length = next(iter(self.adapter_dict[segment_name].values()))
            elif seg_type == longbow.utils.constants.RANDOM_SEGMENT_NAME:
                segment_length = 0
            else:
                # We should never get here:
                raise RuntimeError(f"ERROR: Got an unknown segment type for {segment_name}: {seg_type}")

        # Just one last sanity check:
        if segment_length is None:
            # We should never get here!
            raise RuntimeError(f"ERROR: Got a NoneType segment_length for segment {segment_name}")

        return segment_length

    @staticmethod
    def from_json_file(json_file):
        """Create a LibraryModel instance from the given json file.
        This method will open the file at the given location and use the data in that file to create a LibraryModel."""

        try:
            with open(json_file) as f:
                json_data = json.load(f)
        except FileNotFoundError:
            logger.error(f"File does not exist: {json_file}")
            sys.exit(1)

        return LibraryModel.from_json_obj(json_data)

    @staticmethod
    def from_json_obj(json_data):
        """Create a LibraryModel instance from the given json data.
        This method will use the data in the JSON object to create a LibraryModel."""

        # Get "optional" new model params:
        try:
            named_random_segments = set(json_data["named_random_segments"])
        except KeyError:
            named_random_segments = None

        try:
            coding_region=json_data["coding_region"]
        except KeyError:
            coding_region = None

        try:
            annotation_segments={k: list(v) for k, v in json_data["annotation_segments"].items()}
        except KeyError:
            annotation_segments = None

        m = LibraryModel(
            name=json_data["name"],
            description=json_data["description"],
            version=json_data["version"],
            array_element_structure=tuple(tuple(v) for v in json_data["array_element_structure"]),
            adapters=json_data["adapters"],
            direct_connections={k: set(v) for k, v in json_data["direct_connections"].items()},
            start_element_names=set(json_data["start_element_names"]),
            end_element_names=set(json_data["end_element_names"]),
            # New "optional" model parameters:
            named_random_segments=named_random_segments,
            coding_region=coding_region,
            annotation_segments=annotation_segments,
        )

        return m

    @staticmethod
    def build_pre_configured_model(model_name):
        """Build a pre-configured model based on the given model name.
        If the given name does not appear in self.pre_configured_models, then a KeyError will be thrown."""
        if model_name not in LibraryModel.pre_configured_models:
            raise KeyError(f"Model not found in pre-configured models: {model_name}")

        # New "optional" parameters:
        named_random_segments = None
        if "named_random_segments" in LibraryModel.pre_configured_models[model_name]:
            named_random_segments = LibraryModel.pre_configured_models[model_name]["named_random_segments"]

        coding_region = None
        if "coding_region" in LibraryModel.pre_configured_models[model_name]:
            coding_region = LibraryModel.pre_configured_models[model_name]["coding_region"]

        annotation_segments = None
        if "annotation_segments" in LibraryModel.pre_configured_models[model_name]:
            annotation_segments = LibraryModel.pre_configured_models[model_name]["annotation_segments"]

        return LibraryModel(
            name=model_name,
            description=LibraryModel.pre_configured_models[model_name]["description"],
            version=LibraryModel.pre_configured_models[model_name]["version"],
            array_element_structure=LibraryModel.pre_configured_models[model_name]["array_element_structure"],
            adapters=LibraryModel.pre_configured_models[model_name]["adapters"],
            direct_connections=LibraryModel.pre_configured_models[model_name]["direct_connections"],
            start_element_names=LibraryModel.pre_configured_models[model_name]["start_element_names"],
            end_element_names=LibraryModel.pre_configured_models[model_name]["end_element_names"],
            # New "optional" model parameters:
            named_random_segments=named_random_segments,
            coding_region=coding_region,
            annotation_segments=annotation_segments,
        )

    # TODO: Make an enum for this...
    # Model naming convention:
    #     prefix (mas, isoseq)
    #     modality (bulk, sc, spatial)
    #     input library type (10x5p, 10x3p, slideseq)
    #     umi style (none, single, dual)
    #     plexity (pbkit, internal, none)
    #
    # mas_<num_array_elements>_<modality>_<library_type>_<umi_style>_<plexity>
    #
    # e.g.: mas_15_sc_10x5p_single_none
    pre_configured_models = {
        "mas_15_sc_10x5p_single_none": {
            "description": "The standard MAS-seq 15 array element model.",
            "version": "2.0.1",
            "array_element_structure": (
                # NOTE: the first element doesn't currently have the "A" adapter in this version of the library.
                ("A", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("B", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("C", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("D", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("E", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("F", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("G", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("H", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("I", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("J", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("K", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("L", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("M", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("N", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                # The last element doesn't currently have the "P" adapter in this version of the library:
                ("O", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter", "P"),
            ),
            "adapters": {
                "5p_Adapter": "TCTACACGACGCTCTTCCGATCT",
                # "Poly_A": "A" * 30,
                "Poly_A": {HPR_SEGMENT_TYPE_NAME: ("A", 30)},
                "3p_Adapter": "GTACTCTGCGTTGATACCACTGCTT",
                "SLS": "TTTCTTATATGGG",
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
                "A": {"5p_Adapter"},
                "B": {"5p_Adapter"},
                "C": {"5p_Adapter"},
                "D": {"5p_Adapter"},
                "E": {"5p_Adapter"},
                "F": {"5p_Adapter"},
                "G": {"5p_Adapter"},
                "H": {"5p_Adapter"},
                "I": {"5p_Adapter"},
                "J": {"5p_Adapter"},
                "K": {"5p_Adapter"},
                "L": {"5p_Adapter"},
                "M": {"5p_Adapter"},
                "N": {"5p_Adapter"},
                "O": {"5p_Adapter"},
                "P": {"5p_Adapter"},
                "5p_Adapter": {"CBC"},
                "CBC": {"UMI"},
                "UMI": {"SLS"},
                "SLS": {"cDNA"},
                "cDNA": {"Poly_A"},
            },
            "start_element_names": {"A", "5p_Adapter"},
            "end_element_names": {"Poly_A", "P"},
            "named_random_segments": {"UMI", "cDNA", "CBC"},
            "coding_region": "cDNA",
            "annotation_segments": {
                "UMI": [(longbow.utils.constants.READ_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG),
                        (longbow.utils.constants.READ_RAW_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG)],
                "CBC": [(longbow.utils.constants.READ_BARCODE_TAG, longbow.utils.constants.READ_BARCODE_POS_TAG), (
                longbow.utils.constants.READ_RAW_BARCODE_TAG, longbow.utils.constants.READ_BARCODE_POS_TAG)],
            },
            "deprecated": False,
        },
        "mas_15_sc_10x3p_single_none": {
            "description": "The 3' kit MAS-seq 15 array element model.",
            "version": "2.0.1",
            "array_element_structure": (
                # NOTE: the first element doesn't currently have the "A" adapter in this version of the library.
                ("A", "5p_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("B", "5p_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("C", "5p_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("D", "5p_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("E", "5p_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("F", "5p_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("G", "5p_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("H", "5p_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("I", "5p_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("J", "5p_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("K", "5p_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("L", "5p_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("M", "5p_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("N", "5p_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                # The last element doesn't currently have the "P" adapter in this version of the library:
                ("O", "5p_Adapter", "CBC", "UMI", "Poly_T", "cDNA", "3p_Adapter", "P"),
            ),
            "adapters": {
                "5p_Adapter": "TCTACACGACGCTCTTCCGATCT",
                "3p_Adapter": "CCCATGTACTCTGCGTTGATACCACTGCTT",
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
                "Poly_T": {HPR_SEGMENT_TYPE_NAME: ("T", 30)},
                "CBC": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 16},
                "UMI": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 10},
                "cDNA": RANDOM_SEGMENT_NAME,
            },
            "direct_connections": {
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
                "A": {"5p_Adapter"},
                "B": {"5p_Adapter"},
                "C": {"5p_Adapter"},
                "D": {"5p_Adapter"},
                "E": {"5p_Adapter"},
                "F": {"5p_Adapter"},
                "G": {"5p_Adapter"},
                "H": {"5p_Adapter"},
                "I": {"5p_Adapter"},
                "J": {"5p_Adapter"},
                "K": {"5p_Adapter"},
                "L": {"5p_Adapter"},
                "M": {"5p_Adapter"},
                "N": {"5p_Adapter"},
                "O": {"5p_Adapter"},
                "P": {"5p_Adapter"},
                "5p_Adapter": {"CBC"},
                "CBC": {"UMI"},
                "UMI": {"Poly_T"},
                "Poly_T": {"cDNA"},
                "cDNA": {"3p_Adapter"},
            },
            "start_element_names": {"A", "5p_Adapter"},
            "end_element_names": {"P"},
            "named_random_segments": {"UMI", "cDNA", "CBC"},
            "coding_region": "cDNA",
            "annotation_segments": {
                "UMI": [(longbow.utils.constants.READ_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG),
                        (longbow.utils.constants.READ_RAW_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG)],
                "CBC": [(longbow.utils.constants.READ_BARCODE_TAG, longbow.utils.constants.READ_BARCODE_POS_TAG), (
                longbow.utils.constants.READ_RAW_BARCODE_TAG, longbow.utils.constants.READ_BARCODE_POS_TAG)],
            },
            "deprecated": False,
        },
        "mas_15_bulk_10x5p_single_internal": {
            "description": "A MAS-seq 15 array element model with a 10 base index just before the 3' adapter for bulk "
                           "sequencing.",
            "version": "1.0.1",
            "array_element_structure": (
                ("A", "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter"),
                ("B", "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter"),
                ("C", "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter"),
                ("D", "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter"),
                ("E", "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter"),
                ("F", "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter"),
                ("G", "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter"),
                ("H", "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter"),
                ("I", "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter"),
                ("J", "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter"),
                ("K", "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter"),
                ("L", "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter"),
                ("M", "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter"),
                ("N", "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter"),
                ("O", "5p_Adapter", "UMI", "SLS", "cDNA", "Poly_A", "sample_index", "3p_Adapter", "P"),
            ),
            "adapters": {
                "5p_Adapter": "TCTACACGACGCTCTTCCGATCT",
                "Poly_A": {HPR_SEGMENT_TYPE_NAME: ("A", 30)},
                "3p_Adapter": "CTCTGCGTTGATACCACTGCTT",
                "SLS": "TTTCTTATATGGG",
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
                "sample_index": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 10},
                "UMI": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 10},
                "cDNA": RANDOM_SEGMENT_NAME,
            },
            "direct_connections": {
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
                "A": {"5p_Adapter"},
                "B": {"5p_Adapter"},
                "C": {"5p_Adapter"},
                "D": {"5p_Adapter"},
                "E": {"5p_Adapter"},
                "F": {"5p_Adapter"},
                "G": {"5p_Adapter"},
                "H": {"5p_Adapter"},
                "I": {"5p_Adapter"},
                "J": {"5p_Adapter"},
                "K": {"5p_Adapter"},
                "L": {"5p_Adapter"},
                "M": {"5p_Adapter"},
                "N": {"5p_Adapter"},
                "O": {"5p_Adapter"},
                "P": {"5p_Adapter"},
                "5p_Adapter": {"UMI"},
                "UMI": {"SLS"},
                "SLS": {"cDNA"},
                "cDNA": {"Poly_A"},
                "Poly_A": {"sample_index"},
                "sample_index": {"3p_Adapter"},
            },
            "start_element_names": {"A", "5p_Adapter"},
            "end_element_names": {"Poly_A", "P"},
            "named_random_segments": {"UMI", "cDNA", "sample_index"},
            "coding_region": "cDNA",
            "annotation_segments": {
                "UMI": [(longbow.utils.constants.READ_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG),
                        (longbow.utils.constants.READ_RAW_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG)],
                "sample_index": [(longbow.utils.constants.READ_DEMUX_TAG, longbow.utils.constants.READ_DEMUX_POS_TAG)],
            },
            "deprecated": False,
        },
        "mas_10_sc_10x5p_single_none": {
            "description": "The MAS-seq 10 array element model.",
            "version": "2.0.1",
            "array_element_structure": (
                ("Q", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("C", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("M", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("I", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("O", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("J", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("B", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("D", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                ("K", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter"),
                # The last element may not currently have the "R" adapter in this version of the library:
                ("H", "5p_Adapter", "CBC", "UMI", "SLS", "cDNA", "Poly_A", "3p_Adapter", "R"),
            ),
            "adapters": {
                "5p_Adapter": "TCTACACGACGCTCTTCCGATCT",
                "Poly_A": {HPR_SEGMENT_TYPE_NAME: ("A", 30)},
                "3p_Adapter": "GTACTCTGCGTTGATACCACTGCTT",
                "SLS": "TTTCTTATATGGG",
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
                "B": {"5p_Adapter"},
                "C": {"5p_Adapter"},
                "D": {"5p_Adapter"},
                "H": {"5p_Adapter"},
                "I": {"5p_Adapter"},
                "J": {"5p_Adapter"},
                "K": {"5p_Adapter"},
                "M": {"5p_Adapter"},
                "O": {"5p_Adapter"},
                "Q": {"5p_Adapter"},
                "R": {"5p_Adapter"},
                "5p_Adapter": {"CBC"},
                "CBC": {"UMI"},
                "UMI": {"SLS"},
                "SLS": {"cDNA"},
                "cDNA": {"Poly_A"},
            },
            "start_element_names": {"Q", "5p_Adapter"},
            "end_element_names": {"Poly_A", "R"},
            "named_random_segments": {"UMI", "cDNA", "CBC"},
            "coding_region": "cDNA",
            "annotation_segments": {
                "UMI": [(longbow.utils.constants.READ_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG),
                        (longbow.utils.constants.READ_RAW_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG)],
                "CBC": [(longbow.utils.constants.READ_BARCODE_TAG, longbow.utils.constants.READ_BARCODE_POS_TAG), (
                    longbow.utils.constants.READ_RAW_BARCODE_TAG, longbow.utils.constants.READ_BARCODE_POS_TAG)],
            },
            "deprecated": False,
        },
        "mas_15_spatial_slide-seq_single_none": {
            # The slide-seq model is:
            #
            #                 |-----5p_Adapter---->        |--splitter------>               |------Poly_T---------------->                  |--------5p_Adapter----------|                         # noqa
            # AGCTTACTTGTGAAGACTACACGACGCTCTTCCGATCTNNNNNNNNTCTTCAGCGTTCCCGAGANNNNNNNNNNNNNVVTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVNNNNNNNNNNNNNNNNNCCCATGTACTCTGCGTTGATACCACTGCTTACTTGTAAGCTGTCTA...      # noqa
            # |------A------->                      <------|                  <-----------|                                 <----cDNA-------|                              |-------B------>         # noqa
            #                                          V                           V
            #                                    Spatial Barcode 2         Spatial Barcode 1
            "description": "The Slide-seq 15 array element model.",
            "version": "2.0.1",
            "array_element_structure": (
                ("A", "5p_Adapter", "SBC2", "SLS2", "SBC1", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("B", "5p_Adapter", "SBC2", "SLS2", "SBC1", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("C", "5p_Adapter", "SBC2", "SLS2", "SBC1", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("D", "5p_Adapter", "SBC2", "SLS2", "SBC1", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("E", "5p_Adapter", "SBC2", "SLS2", "SBC1", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("F", "5p_Adapter", "SBC2", "SLS2", "SBC1", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("G", "5p_Adapter", "SBC2", "SLS2", "SBC1", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("H", "5p_Adapter", "SBC2", "SLS2", "SBC1", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("I", "5p_Adapter", "SBC2", "SLS2", "SBC1", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("J", "5p_Adapter", "SBC2", "SLS2", "SBC1", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("K", "5p_Adapter", "SBC2", "SLS2", "SBC1", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("L", "5p_Adapter", "SBC2", "SLS2", "SBC1", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("M", "5p_Adapter", "SBC2", "SLS2", "SBC1", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("N", "5p_Adapter", "SBC2", "SLS2", "SBC1", "UMI", "Poly_T", "cDNA", "3p_Adapter"),
                ("O", "5p_Adapter", "SBC2", "SLS2", "SBC1", "UMI", "Poly_T", "cDNA", "3p_Adapter", "P"),
            ),
            "adapters": {
                "5p_Adapter": "TCTACACGACGCTCTTCCGATCT",
                "3p_Adapter": "CCCATGTACTCTGCGTTGATACCACTGCTT",
                "SLS2": "TCTTCAGCGTTCCCGAGA",
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
                "Poly_T": {HPR_SEGMENT_TYPE_NAME: ("T", 30)},
                "SBC1": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 6},
                "SBC2": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 8},
                # The UMI might be 7, rather than 9 elements long - not clear from the geneious file.
                "UMI": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 9},
                "cDNA": RANDOM_SEGMENT_NAME,
            },
            "direct_connections": {
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
                "A": {"5p_Adapter"},
                "B": {"5p_Adapter"},
                "C": {"5p_Adapter"},
                "D": {"5p_Adapter"},
                "E": {"5p_Adapter"},
                "F": {"5p_Adapter"},
                "G": {"5p_Adapter"},
                "H": {"5p_Adapter"},
                "I": {"5p_Adapter"},
                "J": {"5p_Adapter"},
                "K": {"5p_Adapter"},
                "L": {"5p_Adapter"},
                "M": {"5p_Adapter"},
                "N": {"5p_Adapter"},
                "O": {"5p_Adapter"},
                "5p_Adapter": {"SBC2"},
                "SBC2": {"SLS2"},
                "SLS2": {"SBC1"},
                "SBC1": {"UMI"},
                "UMI": {"Poly_T"},
                "Poly_T": {"cDNA"},
                "cDNA": {"3p_Adapter"},
            },
            "start_element_names": {"A", "5p_Adapter"},
            "end_element_names": {"P", "3p_Adapter"},
            "named_random_segments": {"UMI", "SBC2", "SBC1", "cDNA"},
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
        "mas_15_bulk_teloprimeV2_single_none": {
            "description": "The MAS15 Teloprime V2 indexed array element model.",
            "version": "2.0.1",
            "array_element_structure": (
                ("A", "TPV2_adapter", "cDNA", "Poly_A", "idx", "rev_bind"),
                ("B", "TPV2_adapter", "cDNA", "Poly_A", "idx", "rev_bind"),
                ("C", "TPV2_adapter", "cDNA", "Poly_A", "idx", "rev_bind"),
                ("D", "TPV2_adapter", "cDNA", "Poly_A", "idx", "rev_bind"),
                ("E", "TPV2_adapter", "cDNA", "Poly_A", "idx", "rev_bind"),
                ("F", "TPV2_adapter", "cDNA", "Poly_A", "idx", "rev_bind"),
                ("G", "TPV2_adapter", "cDNA", "Poly_A", "idx", "rev_bind"),
                ("H", "TPV2_adapter", "cDNA", "Poly_A", "idx", "rev_bind"),
                ("I", "TPV2_adapter", "cDNA", "Poly_A", "idx", "rev_bind"),
                ("J", "TPV2_adapter", "cDNA", "Poly_A", "idx", "rev_bind"),
                ("K", "TPV2_adapter", "cDNA", "Poly_A", "idx", "rev_bind"),
                ("L", "TPV2_adapter", "cDNA", "Poly_A", "idx", "rev_bind"),
                ("M", "TPV2_adapter", "cDNA", "Poly_A", "idx", "rev_bind"),
                ("N", "TPV2_adapter", "cDNA", "Poly_A", "idx", "rev_bind"),
                ("O", "TPV2_adapter", "cDNA", "Poly_A", "idx", "rev_bind", "P"),
            ),
            "adapters": {
                "TPV2_adapter": "CTACACGACGCTCTTCCGATCTTGGATTGATATGTAATACGACTCACTATAG",
                "rev_bind": "CTCTGCGTTGATACCACTGCTT",
                "A": "AGCTTACTTGTGAAGAT",
                "B": "ACTTGTAAGCTGTCTAT",
                "C": "ACTCTGTCAGGTCCGAT",
                "D": "ACCTCCTCCTCCAGAAT",
                "E": "AACCGGACACACTTAGT",
                "F": "AGAGTCCAATTCGCAGT",
                "G": "AATCAAGGCTTAACGGT",
                "H": "ATGTTGAATCCTAGCGT",
                "I": "AGTGCGTTGCGAATTGT",
                "J": "AATTGCGTAGTTGGCCT",
                "K": "ACACTTGGTCGCAATCT",
                "L": "AGTAAGCCTTCGTGTCT",
                "M": "ACCTAGATCAGAGCCTT",
                "N": "AGGTATGCCGGTTAAGT",
                "O": "AAGTCACCGGCACCTTT",
                "P": "ATGAAGTGGCTCGAGA",
                "Poly_A": {HPR_SEGMENT_TYPE_NAME: ("A", 30)},
                "idx": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 10},
                "cDNA": RANDOM_SEGMENT_NAME,
            },
            "direct_connections": {
                "A": {"TPV2_adapter"},
                "B": {"TPV2_adapter"},
                "C": {"TPV2_adapter"},
                "D": {"TPV2_adapter"},
                "E": {"TPV2_adapter"},
                "F": {"TPV2_adapter"},
                "G": {"TPV2_adapter"},
                "H": {"TPV2_adapter"},
                "I": {"TPV2_adapter"},
                "J": {"TPV2_adapter"},
                "K": {"TPV2_adapter"},
                "L": {"TPV2_adapter"},
                "M": {"TPV2_adapter"},
                "N": {"TPV2_adapter"},
                "O": {"TPV2_adapter"},
                "TPV2_adapter": {"cDNA"},
                "cDNA": {"Poly_A"},
                "Poly_A": {"idx"},
                "idx": {"rev_bind"},
                "rev_bind": {
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
            },
            "start_element_names": {"A", "TPV2_adapter"},
            "end_element_names": {"P", "rev_bind"},
            "named_random_segments": {"idx", "cDNA"},
            "coding_region": "cDNA",
            "annotation_segments": {
                "idx": [(longbow.utils.constants.READ_INDEX_TAG, longbow.utils.constants.READ_BARCODE_POS_TAG)],
            },
            "deprecated": False,
        },
        "isoseq_1_sc_10x5p_single_none": {
            "description": "Single-cell RNA (without MAS-seq prep).",
            "version": "1.0.1",
            "array_element_structure": (
                ("V", "CBC", "UMI", "B", "cDNA", "Poly_A", "M"),
            ),
            "adapters": {
                "V": "TCTACACGACGCTCTTCCGATCT",
                "Poly_A": {HPR_SEGMENT_TYPE_NAME: ("A", 30)},
                "M": "GTACTCTGCGTTGATACCACTGCTT",
                "B": "TTTCTTATATGGG",
                "CBC": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 16},
                "UMI": {FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME: 10},
                "cDNA": RANDOM_SEGMENT_NAME,
            },
            "direct_connections": {
                "V": {"CBC"},
                "CBC": {"UMI"},
                "UMI": {"B"},
                "B": {"cDNA"},
                "cDNA": {"Poly_A"},
                "Poly_A": {"M"},
            },
            "start_element_names": {"V"},
            "end_element_names": {"M"},
            "named_random_segments": {"UMI", "cDNA", "CBC"},
            "coding_region": "cDNA",
            "annotation_segments": {
                "UMI": [(longbow.utils.constants.READ_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG),
                        (longbow.utils.constants.READ_RAW_UMI_TAG, longbow.utils.constants.READ_UMI_POS_TAG)],
                "CBC": [(longbow.utils.constants.READ_BARCODE_TAG, longbow.utils.constants.READ_BARCODE_POS_TAG), (
                longbow.utils.constants.READ_RAW_BARCODE_TAG, longbow.utils.constants.READ_BARCODE_POS_TAG)],
            },
            "deprecated": False,
        },
    }

