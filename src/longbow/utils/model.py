import sys
import re
import logging

import click_log

from pomegranate import *

import longbow.utils.model_utils
from .model_utils import ModelBuilder

import longbow.utils.constants
from .constants import FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME, HPR_SEGMENT_TYPE_NAME, START_STATE_INDICATOR, END_STATE_INDICATOR

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger(__name__)
click_log.basic_config(logger)


class LibraryModel:
    """Model describing a given library preparation.
    The model can annotate the known sections of a read from the library it describes."""

    def __init__(self,
                 array_model,
                 cdna_model,
                 model_name='hierarchical_model',
                 do_build=True):

        self.name = model_name
        self.description = f"{array_model['description']}, {cdna_model['description']}"
        self.array_version = array_model['version']
        self.cdna_version = array_model['version']

        self.array_model = array_model
        self.cdna_model = cdna_model

        self.adapter_dict = {**array_model['adapters'], **cdna_model['adapters']}

        self.named_random_segments = cdna_model['named_random_segments']
        self.coding_region = cdna_model['coding_region']
        self.annotation_segments = cdna_model['annotation_segments']

        self.hmm = None
        self.key_adapters = self.array_model['structure']
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
        if self.cdna_model['annotation_segments']:
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
        return len(self.array_model['structure'])

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
            if n in self.array_model['structure']:

                # If this is our first segment, we should allow for the possibility that our array begins
                # somewhere after the first element.  We must find the starting point:
                if not found_first_key:
                    while n != self.array_model['structure'][key_adapter_indx]:
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
                if key_adapter_indx < len(self.array_model['structure']) and (n == self.array_model['structure'][key_adapter_indx]):
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
                (self.array_model['structure'][first_key_adapter_index] == self.array_model['structure'][-1]):
            is_valid = False

        return is_valid, num_key_adapters_found, first_key_adapter_index

    def get_all_node_names(self):
        """Return a list containing the names of all nodes in this model."""

        # Let's be exhaustive here just to make sure we don't miss anything:
        all_node_names = set()

        for e in self.array_model['structure']:
            all_node_names.add(e)

        for e in self.cdna_model['structure']:
            all_node_names.add(e)

        return sorted(list(all_node_names))

    def build(self):
        # Initiate model
        random_states, random_model = self._create_random_repeat_model()
        self.hmm = random_model

        # Add cDNA model
        cdna_states, cdna_hmm = self._create_cdna_model()
        self.hmm.add_model(cdna_hmm)

        # Add array model
        array_states, array_hmm = self._create_array_model()
        self.hmm.add_model(array_hmm)

        # Make a dict of all available states for easy lookup and connection
        all_states = {**random_states, **cdna_states, **array_states}

        # Connect to model start and stop states
        terminal_states = self.array_model['structure'] + self.cdna_model['structure']
        for terminal_state in terminal_states:
            self.hmm.add_transition(self.hmm.start, all_states[f'{terminal_state}-start'], 1.0/len(terminal_states))
            self.hmm.add_transition(all_states[f'{terminal_state}-end'], self.hmm.end, 1.0/len(terminal_states))

        # Connect array end states to random model
        for i in range(len(self.array_model['structure'])-1):
            self.hmm.add_transition(
                all_states[f'{self.array_model["structure"][i]}-end'],
                all_states['random:RDA'],
                1.0
            )

        # Connect random model to array start states
        for i in range(1, len(self.array_model['structure'])):
            self.hmm.add_transition(
                all_states['random:RDB'],
                all_states[f'{self.array_model["structure"][i]}-start'],
                1.0/(len(self.array_model['structure'])-1)
            )

        self.hmm.bake()

    def _create_random_repeat_model(self):
        model = ModelBuilder.make_random_repeat_model()

        states = dict()
        for s in model.states:
            states[s.name] = s

        return states, model

    def _create_array_model(self):
        model = None

        for adapter_name in self.array_model['structure']:
            adapter_hmm = ModelBuilder.make_global_alignment_model(self.array_model['adapters'][adapter_name], adapter_name)

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
        model = HiddenMarkovModel(name='cdna-model')

        # Make the individual adapter models
        for adapter_name, adapter_def in self.cdna_model['adapters'].items():
            adapter_hmm = None

            if type(adapter_def) is str:
                if adapter_name in self.cdna_model['named_random_segments']:
                    adapter_hmm = ModelBuilder.make_named_random_model(adapter_name)
                else:
                    adapter_hmm = ModelBuilder.make_global_alignment_model(self.cdna_model['adapters'][adapter_name], adapter_name)

            elif type(adapter_def) is dict:
                segment_type = list(adapter_def.keys())[0]

                if segment_type == FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME:
                    adapter_hmm = ModelBuilder.make_fixed_length_random_segment(adapter_name, list(adapter_def.values())[0])
                elif segment_type == HPR_SEGMENT_TYPE_NAME:
                    base, hpr_length = list(adapter_def.values())[0]
                    adapter_hmm = ModelBuilder.make_homopolymer_repeat_model(adapter_name, base, hpr_length)
                else:
                    raise RuntimeError(f"Unknown special model type: {segment_type}")

            if model is None:
                model = adapter_hmm
            else:
                model.add_model(adapter_hmm)

        # Intermediate bake to facilitate searching for states in the model.
        model.bake(merge="None")

        # Connect components
        for i in range(0, len(self.cdna_model['structure'])-1):
            for j in range(i+1, len(self.cdna_model['structure'])):
                adapter_name_i = self.cdna_model['structure'][i]
                adapter_name_j = self.cdna_model['structure'][j]

                prob = 0.9 if j - i == 1 else 0.1 / (len(self.cdna_model['structure']) - i - 2)
                ModelBuilder.connect_terminals(model, adapter_name_i, adapter_name_j, prob)

        cstart = ModelBuilder.find_state(model, 'cdna-model-start')
        sstart = ModelBuilder.find_state(model, f'{self.cdna_model["structure"][0]}-start')
        send = ModelBuilder.find_state(model, f'{self.cdna_model["structure"][-1]}-end')
        cend = ModelBuilder.find_state(model, 'cdna-model-end')

        model.add_transition(cstart, sstart, 1.0)
        model.add_transition(send, cend, 1.0)

        model.bake(merge="None")

        states = dict()
        for s in model.states:
            states[s.name] = s

        return states, model

    # TODO: fix this
    def to_json(self, outfile=None, indent=4):
        """Serialize this model to a json object and return that json object.
        If outfile is not none, will write the json object to the given file path."""

        model_data = {
            "name": self.name,
            "description": self.description,
            # "version": self.version,
            # "array_element_structure": self.array_element_structure,
            # "adapters": self.adapter_dict,
            # "direct_connections": {k: list(v) for k, v in self.direct_connections_dict.items()},
            # "start_element_names": list(self.start_element_names),
            # "end_element_names": list(self.end_element_names),
        }

        if outfile:
            with open(outfile, 'w') as f:
                json.dump(model_data, f, indent=indent)

        return json.dumps(model_data, indent=indent)

    @staticmethod
    def has_prebuilt_model(model_name):
        (array_model_name, cdna_model_name) = re.split('\+', model_name, 2)

        if array_model_name not in ModelBuilder.pre_configured_array_models.keys():
            return False

        if cdna_model_name not in ModelBuilder.pre_configured_cdna_models.keys():
            return False

        return True

    @staticmethod
    def build_pre_configured_model(model_name):
        (array_model_name, cdna_model_name) = re.split('\+', model_name, 2)

        lb = LibraryModel(
            array_model=ModelBuilder.pre_configured_array_models[array_model_name], 
            cdna_model=ModelBuilder.pre_configured_cdna_models[cdna_model_name],
            model_name=model_name
        )

        return lb