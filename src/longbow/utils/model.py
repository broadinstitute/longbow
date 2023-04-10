import json
import logging
import re
import sys

import longbow.utils.constants
import longbow.utils.model_utils

from ..meta import VERSION
from .constants import (
    END_STATE_INDICATOR,
    FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME,
    HPR_SEGMENT_TYPE_NAME,
    MODEL_DESC_DELIMITER,
    RANDOM_SEGMENT_NAME,
    START_STATE_INDICATOR,
)
from .model_utils import ModelBuilder

logger = logging.getLogger(__name__)

# NOTE: The backslash here is because the delimiter itself is a plus sign.
MODEL_NAME_REGEX = re.compile(f"\\{MODEL_DESC_DELIMITER}")


class LibraryModel:
    """Model describing a given library preparation.
    The model can annotate the known sections of a read from the library it describes.
    """

    def __init__(
        self,
        array_model,
        cdna_model,
        model_name=None,
        model_description=None,
        do_build=True,
    ):
        if not model_name:
            self.name = array_model["name"] + MODEL_DESC_DELIMITER + cdna_model["name"]
        else:
            self.name = model_name

        if not model_description:
            self.description = (
                f"{array_model['description']}, {cdna_model['description']}"
            )
        else:
            self.description = model_description

        self.array_version = array_model["version"]
        self.cdna_version = cdna_model["version"]

        self.array_model = array_model
        self.cdna_model = cdna_model

        self.adapter_dict = {**array_model["adapters"], **cdna_model["adapters"]}

        self.named_random_segments = set(cdna_model["named_random_segments"])
        self.coding_region = cdna_model["coding_region"]
        self.annotation_segments = cdna_model["annotation_segments"]

        self.hmm = None
        self.key_adapters = self.array_model["structure"]
        self.key_adapter_set = set(self.key_adapters)

        # Reconstruct array element structure in a way that's compatible with the segment tool
        l1 = []
        for i in range(len(self.array_model["structure"]) - 1):
            ae = self.array_model["structure"][i]
            l1.append([ae] + self.cdna_model["structure"])

        ae = self.array_model["structure"][-1]
        l1[-1].append(ae)

        l2 = []
        for x in l1:
            l2.append(tuple(x))

        self.array_element_structure = tuple(l2)

        if do_build:
            self.build()

    @property
    def has_umi_annotation(self):
        return self._has_annotation(
            longbow.utils.constants.READ_RAW_UMI_TAG
        ) or self._has_annotation(longbow.utils.constants.READ_UMI_TAG)

    @property
    def has_cell_barcode_annotation(self):
        return self._has_annotation(
            longbow.utils.constants.READ_RAW_BARCODE_TAG
        ) or self._has_annotation(longbow.utils.constants.READ_BARCODE_TAG)

    def _has_annotation(self, annotation_tag):
        if self.cdna_model["annotation_segments"]:
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
        return len(self.array_model["structure"])

    def annotate(self, seq):
        """Annotate the given segment using this model."""
        logp, path = self.hmm.viterbi(seq)

        apath = []
        ppath = []

        cigar = []
        cur_adapter_name = ""
        cur_op = ""
        cur_op_len = 0

        for idx, state in path:
            apath.append(re.sub(r"\d+$", "", state.name))

            if not state.name.endswith(
                START_STATE_INDICATOR
            ) and not state.name.endswith(END_STATE_INDICATOR):
                adapter_name, op = re.split(r"[:-]", re.sub(r"\d+$", "", state.name), 2)

                if adapter_name != cur_adapter_name:
                    if cur_adapter_name != "":
                        cigar.append(f"{cur_op}{cur_op_len}")
                        ppath.append(f'{cur_adapter_name}:{"".join(cigar)}')

                    cigar = []
                    cur_adapter_name = adapter_name
                    cur_op = op
                    cur_op_len = 0

                if op != cur_op:
                    if cur_op != "":
                        cigar.append(f"{cur_op}{cur_op_len}")

                    cur_op = op
                    cur_op_len = 0

                cur_op_len += 1

        cigar.append(f"{cur_op}{cur_op_len}")
        ppath.append(f'{cur_adapter_name}:{"".join(cigar)}')

        return logp, ppath

    def validate_segment_order(self, segments, allow_missing_first_adapter=True):
        """Validate the order of the given segments against the expected order in this model.

        Returns: (True|False, # key adapters found, first key adapter index)"""

        # Iterate through our given segment names and check if they occur in order:
        ordered_segment_names = tuple([s.name for s in segments])

        num_key_adapters_found = 0
        key_adapter_indx = 0
        first_key_adapter_index = 0
        found_first_key = False
        for n in ordered_segment_names:
            # Ignore all segment names that do not characterize our library:
            if n in self.array_model["structure"]:
                # If this is our first segment, we should allow for the possibility that our array begins
                # somewhere after the first element.  We must find the starting point:
                if not found_first_key:
                    while n != self.array_model["structure"][key_adapter_indx]:
                        key_adapter_indx += 1
                    first_key_adapter_index = key_adapter_indx
                    found_first_key = True

                    # TODO: This can be eliminated for newer datasets, but is here because we're still testing with
                    #       older data that does not have the first MAS-seq overhang.  The model has already been
                    #       updated to handle these segments.
                    # TODO: Replace all references to "5p_Adapter" with references to the models themselves.
                    if (
                        allow_missing_first_adapter
                        and (key_adapter_indx == 1)
                        and (ordered_segment_names[0] == "5p_Adapter")
                    ):
                        num_key_adapters_found += 1
                        first_key_adapter_index = 0

                # Check our key segments here:
                # NOTE: it's possible that we can start matching in the middle of an array, but at the start of a read,
                #       and so we have to do a bounds check here:
                if key_adapter_indx < len(self.array_model["structure"]) and (
                    n == self.array_model["structure"][key_adapter_indx]
                ):
                    key_adapter_indx += 1
                    num_key_adapters_found += 1
                else:
                    # This read does not conform to the model!
                    return False, num_key_adapters_found, first_key_adapter_index

        # If we've made it here and we have seen at least 1 key adapter, then we have a valid array:
        is_valid = True if num_key_adapters_found > 0 else False

        # Except!  If we have an array singleton and it's the terminal overhang adapter,
        # then we do NOT have a valid array:
        if (
            is_valid
            and (num_key_adapters_found == 1)
            and (
                self.array_model["structure"][first_key_adapter_index]
                == self.array_model["structure"][-1]
            )
        ):
            is_valid = False

        return is_valid, num_key_adapters_found, first_key_adapter_index

    def get_all_node_names(self):
        """Return a list containing the names of all nodes in this model."""

        # Let's be exhaustive here just to make sure we don't miss anything:
        all_node_names = set()

        for e in self.array_model["structure"]:
            all_node_names.add(e)

        for e in self.cdna_model["structure"]:
            all_node_names.add(e)

        return sorted(list(all_node_names))

    @staticmethod
    def _warn_if_deprecated(m, model_type):
        try:
            if m["deprecated"]:
                logger.warning(f"{model_type} Model {m.name} is deprecated.")
        except KeyError:
            logger.warning(f"{model_type} Model is missing the 'deprecated' field.")
            pass

    def build(self):
        # Validate our model here so we can go through and just worry about creating it:
        self.validate_model()

        # Warn if either primary sub-model is deprecated:
        LibraryModel._warn_if_deprecated(self.array_model, "Array")
        LibraryModel._warn_if_deprecated(self.cdna_model, "cDNA")

        # Initiate model
        random_states, random_model = LibraryModel._create_random_repeat_model()
        self.hmm = random_model

        # Add cDNA model
        cdna_states, cdna_hmm = self._create_cdna_model()
        self.hmm.add_model(cdna_hmm)

        # Add array model
        array_states, array_hmm = self._create_array_model()
        self.hmm.add_model(array_hmm)

        # Make a dict of all available states for easy lookup and connection
        all_states = {**random_states, **cdna_states, **array_states}

        # Connect random start to array start and cDNA start
        self.hmm.add_transition(
            all_states["random-start"],
            all_states[f"{self.cdna_model['structure'][0]}-start"],
            0.5,
        )
        self.hmm.add_transition(
            all_states["random-start"],
            all_states[f"{self.array_model['structure'][0]}-start"],
            0.5,
        )

        # Connect random:RDA to appropriate array starts
        for array_state in self.array_model["structure"][
            1 : len(self.array_model["structure"])
        ]:
            self.hmm.add_transition(
                all_states["random:RDA"], all_states[f"{array_state}-start"], 0.05
            )

        # Connect random:RDA to non-random cDNA starts
        for cdna_state in self.cdna_model["structure"][
            1 : len(self.cdna_model["structure"])
        ]:
            if cdna_state not in self.named_random_segments:
                self.hmm.add_transition(
                    all_states["random:RDA"], all_states[f"{cdna_state}-start"], 0.05
                )

        # Connect all array ends to random end and cDNA start:
        for array_state in self.array_model["structure"]:
            length = len(self.array_model["adapters"][array_state])

            for op in ["M", "I", "D"]:
                self.hmm.add_transition(
                    all_states[f"{array_state}:{op}{length}"],
                    all_states["random:RDB"],
                    0.01,
                )
                self.hmm.add_transition(
                    all_states[f"{array_state}:{op}{length}"],
                    all_states["random-end"],
                    0.01,
                )

                # NOTE: the last array adapter does not get to have a connection to the cDNA start:
                if array_state != self.array_model["structure"][-1]:
                    self.hmm.add_transition(
                        all_states[f"{array_state}:{op}{length}"],
                        all_states[f'{self.cdna_model["structure"][0]}-start'],
                        1.0,
                    )

        # Connect all cDNA ends to array starts
        for array_state in self.array_model["structure"]:
            cdna_state = self.cdna_model["structure"][-1]
            length = len(self.cdna_model["adapters"][cdna_state])

            for op in ["M", "I", "D"]:
                self.hmm.add_transition(
                    all_states[f"{cdna_state}:{op}{length}"],
                    all_states["random:RDB"],
                    0.01,
                )
                self.hmm.add_transition(
                    all_states[f"{cdna_state}:{op}{length}"],
                    all_states[f"{array_state}-start"],
                    1.0 / (len(self.array_model["structure"])),
                )

        # Connect the cDNA adapter ends to random end
        for i in range(len(self.cdna_model["structure"])):
            adapter_name = self.cdna_model["structure"][i]
            adapter_def = self.cdna_model["adapters"][adapter_name]

            next_adapter_name = (
                None
                if i == len(self.cdna_model["structure"]) - 1
                else self.cdna_model["structure"][i + 1]
            )

            if type(adapter_def) is str:
                if adapter_def == RANDOM_SEGMENT_NAME:
                    self.hmm.add_transition(
                        all_states[f"{adapter_name}-end"],
                        all_states["random:RDB"],
                        0.01,
                    )
                    self.hmm.add_transition(
                        all_states[f"{adapter_name}-end"],
                        all_states["random-end"],
                        0.01,
                    )

                    if next_adapter_name is not None:
                        self.hmm.add_transition(
                            all_states[f"{adapter_name}-end"],
                            all_states[f"{next_adapter_name}-start"],
                            1.0,
                        )
                else:
                    for op in ["M", "I", "D"]:
                        length = len(self.cdna_model["adapters"][adapter_name])
                        self.hmm.add_transition(
                            all_states[f"{adapter_name}:{op}{length}"],
                            all_states["random:RDB"],
                            0.01,
                        )
                        self.hmm.add_transition(
                            all_states[f"{adapter_name}:{op}{length}"],
                            all_states["random-end"],
                            0.01,
                        )

                        if next_adapter_name is not None:
                            self.hmm.add_transition(
                                all_states[f"{adapter_name}:{op}{length}"],
                                all_states[f"{next_adapter_name}-start"],
                                1.0,
                            )

            elif type(adapter_def) is dict:
                segment_type = list(adapter_def.keys())[0]

                if segment_type == FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME:
                    length = list(adapter_def.values())[0]
                    self.hmm.add_transition(
                        all_states[f"{adapter_name}:M{length}"],
                        all_states["random:RDB"],
                        0.01,
                    )
                    self.hmm.add_transition(
                        all_states[f"{adapter_name}:M{length}"],
                        all_states["random-end"],
                        0.01,
                    )

                    if next_adapter_name is not None:
                        self.hmm.add_transition(
                            all_states[f"{adapter_name}:M{length}"],
                            all_states[f"{next_adapter_name}-start"],
                            1.0,
                        )
                elif segment_type == HPR_SEGMENT_TYPE_NAME:
                    base, hpr_length = list(adapter_def.values())[0]
                    for op in ["M", "I"]:
                        self.hmm.add_transition(
                            all_states[f"{adapter_name}:{op}{hpr_length}"],
                            all_states["random:RDB"],
                            0.01,
                        )
                        self.hmm.add_transition(
                            all_states[f"{adapter_name}:{op}{hpr_length}"],
                            all_states["random-end"],
                            0.01,
                        )

                        if next_adapter_name is not None:
                            self.hmm.add_transition(
                                all_states[f"{adapter_name}:{op}{hpr_length}"],
                                all_states[f"{next_adapter_name}-start"],
                                1.0,
                            )

        self.hmm.bake(merge="None")

    @staticmethod
    def _create_random_repeat_model():
        model = ModelBuilder.make_random_repeat_model()

        states = dict()
        for s in model.states:
            states[s.name] = s

        return states, model

    def _create_array_model(self):
        model = None

        for adapter_name in self.array_model["structure"]:
            adapter_hmm = ModelBuilder.make_global_alignment_model(
                self.array_model["adapters"][adapter_name], adapter_name
            )

            if model is None:
                model = adapter_hmm
            else:
                model.add_model(adapter_hmm)

        model.bake(merge="None")

        states = dict()
        for s in model.states:
            states[s.name] = s

        return states, model

    def _create_cdna_model(self):
        model = None

        # Make the individual adapter models
        for adapter_name, adapter_def in self.cdna_model["adapters"].items():
            adapter_hmm = None

            if type(adapter_def) is str:
                if adapter_name in self.named_random_segments:
                    adapter_hmm = ModelBuilder.make_named_random_model(adapter_name)
                else:
                    adapter_hmm = ModelBuilder.make_global_alignment_model(
                        self.cdna_model["adapters"][adapter_name], adapter_name
                    )

            elif type(adapter_def) is dict:
                segment_type = list(adapter_def.keys())[0]

                if segment_type == FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME:
                    adapter_hmm = ModelBuilder.make_fixed_length_random_segment(
                        adapter_name, list(adapter_def.values())[0]
                    )
                elif segment_type == HPR_SEGMENT_TYPE_NAME:
                    base, hpr_length = list(adapter_def.values())[0]
                    adapter_hmm = ModelBuilder.make_homopolymer_repeat_model(
                        adapter_name, base, hpr_length
                    )
                else:
                    raise RuntimeError(f"Unknown special model type: {segment_type}")

            if model is None:
                model = adapter_hmm
            else:
                model.add_model(adapter_hmm)

        # Intermediate bake to facilitate searching for states in the model.
        model.bake(merge="None")

        states = dict()
        for s in model.states:
            states[s.name] = s

        # Create connections:
        # For "soft" adapters (where sequence is unknown), just connect to the next thing
        # "hard" adapters transitions are defined elsewhere
        for i in range(0, len(self.cdna_model["structure"]) - 1):
            adapter_name_i = self.cdna_model["structure"][i]
            adapter_name_j = self.cdna_model["structure"][i + 1]

            adapter_def = self.cdna_model["adapters"][adapter_name_i]

            if type(adapter_def) is dict:
                segment_type = list(adapter_def.keys())[0]

                if segment_type == FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME:
                    length = list(adapter_def.values())[0]
                    model.add_transition(
                        states[f"{adapter_name_i}:M{length}"],
                        states[f"{adapter_name_j}-start"],
                        1.0,
                    )

        model.bake(merge="None")

        return states, model

    @staticmethod
    def _assert_valid_model_json_data(json_data):
        required_fields = ["array", "cdna", "name", "description"]
        missing_fields = []
        for f in required_fields:
            if f not in json_data:
                missing_fields.append(f)
        if len(missing_fields) > 0:
            message = (
                f"ERROR: JSON model is missing one or more required fields: {','.join(missing_fields)}.  "
                f"As of Longbow v{VERSION} all JSON model definitions must have the following fields: "
                f"{','.join(required_fields)}"
            )
            logger.critical(message)
            raise RuntimeError(message)

        if len(json_data) > len(required_fields):
            logger.warning(
                "Model JSON data has additional fields that are not required.  "
                "These fields will be ignored: "
            )

    @staticmethod
    def from_json_file(json_file):
        """Create a LibraryModel instance from the given json file.
        This method will open the file at the given location and use the data in that file to create a LibraryModel.
        """

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

        LibraryModel._assert_valid_model_json_data(json_data)

        return LibraryModel(
            array_model=json_data["array"],
            cdna_model=json_data["cdna"],
            model_name=json_data["name"],
            model_description=json_data["description"],
        )

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
            elif (
                seg_type
                == longbow.utils.constants.FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME
            ):
                segment_length = next(iter(self.adapter_dict[segment_name].values()))
            elif seg_type == longbow.utils.constants.RANDOM_SEGMENT_NAME:
                segment_length = 0
            else:
                # We should never get here:
                raise RuntimeError(
                    f"ERROR: Got an unknown segment type for {segment_name}: {seg_type}"
                )

        # Just one last sanity check:
        if segment_length is None:
            # We should never get here!
            raise RuntimeError(
                f"ERROR: Got a NoneType segment_length for segment {segment_name}"
            )

        return segment_length

    def dump_as_dotfile_simple(self, out_file_name=None):
        """Dump this LibraryModel to disk as a simple dot file for visualization."""

        if out_file_name is None:
            out_file_name = f"longbow_model_{self.name}.Av{self.array_version}.Cv{self.cdna_version}.simple.dot"

        # Get the model as a dictionary so we can handle it
        model_dict = self.hmm.to_dict()

        # Save a map here of node names so we can make some subgraphs:
        named_node_state_maps = {n: [] for n in self.get_all_node_names()}
        named_node_state_maps[RANDOM_SEGMENT_NAME] = []

        with open(out_file_name, "w") as f:
            f.write(f"digraph {self.name} {{\n\n")

            # Write out the nodes:
            f.write("    // States:\n")
            states = set()
            for node in model_dict["states"]:
                name = LibraryModel._get_state_base_name(node["name"])
                states.add(name)
            for state in states:
                f.write(f'    "{state}";\n')
            f.write("\n")

            # Write out the edges:
            f.write("    // Transitions:\n")
            transitions = dict()
            for edge in model_dict["edges"]:
                source = model_dict["states"][edge[0]]
                dest = model_dict["states"][edge[1]]
                weight = edge[2]

                source_name = LibraryModel._get_state_base_name(source["name"])
                dest_name = LibraryModel._get_state_base_name(dest["name"])

                # TODO: This may be too aggressive a filtering criterion.
                if source_name == dest_name:
                    continue

                transition_descriptor = f'"{source_name}" -> "{dest_name}"'

                transitions[transition_descriptor] = weight
            for t, w in transitions.items():
                f.write(f"    {t} [weight={w:.3f}];\n")

            f.write("}\n")

    def dump_as_dotfile(self, out_file_name=None, do_subgraphs=True):
        """Dump this LibraryModel to disk as a dot file for visualization."""

        if out_file_name is None:
            out_file_name = f"longbow_model_{self.name}.Av{self.array_version}.Cv{self.cdna_version}.dot"

        # Get the model as a dictionary so we can handle it
        model_dict = self.hmm.to_dict()

        # Save a map here of node names so we can make some subgraphs:
        named_node_state_maps = {n: [] for n in self.get_all_node_names()}
        named_node_state_maps[RANDOM_SEGMENT_NAME] = []

        with open(out_file_name, "w") as f:
            f.write(f"digraph {self.name} {{\n\n")

            # Write out the nodes:
            f.write("    // States:\n")
            if do_subgraphs:
                # Collect our subgraphs:
                for node in model_dict["states"]:
                    name = LibraryModel._get_state_base_name(node["name"])
                    named_node_state_maps[name].append(node["name"])

                # Write out the subgraphs we have found:
                for subgraph_name, node_names in named_node_state_maps.items():
                    f.write(f'    subgraph "{subgraph_name}"')
                    f.write(" {\n")
                    for node_name in node_names:
                        f.write(f'        "{node_name}";\n')
                    f.write("    }\n")
            else:
                for node in model_dict["states"]:
                    f.write(f"    \"{node['name']}\";\n")

            f.write("\n")

            # Write out the edges:
            f.write("    // Transitions:\n")
            for edge in model_dict["edges"]:
                source = model_dict["states"][edge[0]]
                dest = model_dict["states"][edge[1]]
                weight = edge[2]

                f.write(
                    f"    \"{source['name']}\" -> \"{dest['name']}\" [weight={weight:.3f}];\n"
                )

            f.write("}\n")

    @staticmethod
    def _assert_required_model_fields_are_present(m, required_fields, model_type):
        missing_fields = []
        for f in required_fields:
            if f not in m:
                missing_fields.append(f)
        if len(missing_fields) > 0:
            message = (
                f"ERROR: {model_type} model does not have the following required fields: "
                f"{','.join(missing_fields)}.  All {model_type} models must have the following fields defined: "
                f"{','.join(required_fields)}.  "
                f"Fields may be defined as `None` if they are not logically required."
            )
            logger.critical(message)
            raise RuntimeError(message)

    @staticmethod
    def _assert_model_segments_are_in_collection(
        model, segments_name, collection_name, model_type
    ):
        not_present = []
        if type(model[segments_name]) is str:
            if model[segments_name] not in model[collection_name]:
                not_present.append(model[segments_name])
        else:
            for s in model[segments_name]:
                if s not in model[collection_name]:
                    not_present.append(s)

        if len(not_present) > 0:
            message = (
                f"ERROR: {model_type} model has the following segments defined in {segments_name}, "
                f"but they are not present in {collection_name}: "
                f"{','.join(not_present)}."
            )

            logger.critical(message)
            raise RuntimeError(message)

    def validate_model(self):
        """Ensure that the configuration of this model is well-formed and conforms to the semantics that we require."""

        ##############################################################################
        # Checks on submodel definition structure:
        required_array_model_fields = [
            "description",
            "version",
            "structure",
            "adapters",
            "deprecated",
            "name",
        ]
        required_cdna_model_fields = required_array_model_fields + [
            "named_random_segments",
            "coding_region",
            "annotation_segments",
        ]
        # Make sure our models have required fields:
        LibraryModel._assert_required_model_fields_are_present(
            self.array_model, required_array_model_fields, "Array"
        )
        LibraryModel._assert_required_model_fields_are_present(
            self.cdna_model, required_cdna_model_fields, "cDNA"
        )

        ##############################################################################
        # Checks on general sub-model structure:

        # Array model:
        LibraryModel._assert_model_segments_are_in_collection(
            self.array_model, "adapters", "structure", "Array"
        )
        LibraryModel._assert_model_segments_are_in_collection(
            self.array_model, "structure", "adapters", "Array"
        )

        # cDNA Model:
        LibraryModel._assert_model_segments_are_in_collection(
            self.cdna_model, "adapters", "structure", "cDNA"
        )
        LibraryModel._assert_model_segments_are_in_collection(
            self.cdna_model, "structure", "adapters", "cDNA"
        )

        for seg_def in [
            "named_random_segments",
            "coding_region",
            "annotation_segments",
        ]:
            LibraryModel._assert_model_segments_are_in_collection(
                self.cdna_model, seg_def, "structure", "cDNA"
            )
            LibraryModel._assert_model_segments_are_in_collection(
                self.cdna_model, seg_def, "adapters", "cDNA"
            )

    def to_json(self, outfile=None, indent=4):
        """Serialize this model to a json object and return that json object.
        If outfile is not none, will write the json object to the given file path."""

        model_data = {
            "name": self.name,
            "description": self.description,
            "array": self.array_model,
            "cdna": self.cdna_model,
        }

        if outfile:
            with open(outfile, "w") as f:
                json.dump(model_data, f, indent=indent)

        return json.dumps(model_data, indent=indent)

    @staticmethod
    def has_prebuilt_model(model_name):
        model_name_pieces = MODEL_NAME_REGEX.split(model_name)

        if len(model_name_pieces) == 2:
            array_model_name, cdna_model_name = model_name_pieces

            if (
                array_model_name
                not in ModelBuilder.pre_configured_models["array"].keys()
            ):
                return False

            if cdna_model_name not in ModelBuilder.pre_configured_models["cdna"].keys():
                return False

            return True

        return False

    @staticmethod
    def build_pre_configured_model(model_name):
        (array_model_name, cdna_model_name) = MODEL_NAME_REGEX.split(model_name, 2)

        lb = LibraryModel(
            array_model=ModelBuilder.pre_configured_models["array"][array_model_name],
            cdna_model=ModelBuilder.pre_configured_models["cdna"][cdna_model_name],
            model_name=model_name,
        )

        return lb

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
