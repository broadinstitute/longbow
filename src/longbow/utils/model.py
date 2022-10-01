import sys
import re
import logging

import click_log

from pomegranate import *

import longbow.utils.model_utils
from .model_utils import ModelBuilder

import longbow.utils.constants
from .constants import FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME, HPR_SEGMENT_TYPE_NAME

logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger(__name__)
click_log.basic_config(logger)


class LibraryModel:
    """Model describing a given library preparation.
    The model can annotate the known sections of a read from the library it describes."""

    def __init__(self,
                 array_model,
                 cdna_model,
                 do_build=True):

        self.array_model = array_model
        self.cdna_model = cdna_model
        self.hmm = None

        if do_build:
            self.build()

    def build(self):
        # Initiate cDNA model
        self.hmm = self._create_cdna_model()

        # Add array model
        self.hmm.add_model(self._create_array_model())
        self.hmm.bake(merge="None")

        for i in range(len(self.array_model['structure'])-1):
            ModelBuilder.connect_terminals(self.hmm, self.array_model['structure'][i], 'cdna-model', 0.1)
            ModelBuilder.connect_terminals(self.hmm, 'cdna-model', self.array_model['structure'][i+1], 0.1)

    def _create_array_model(self):
        model = None

        for adapter_name in self.array_model['structure']:
            adapter_hmm = ModelBuilder.make_global_alignment_model(self.array_model['adapters'][adapter_name], adapter_name)

            if model is None:
                model = adapter_hmm
            else:
                model.add_model(adapter_hmm)

        model.bake(merge="None")
        return model

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

        return model


