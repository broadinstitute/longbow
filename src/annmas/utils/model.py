from pomegranate.callbacks import History, ModelCheckpoint
from pomegranate import *

import matplotlib.pyplot as plt
from matplotlib import transforms

import pysam
import pandas as pd
import numpy as np

import re
import math


def make_global_alignment_model(target, name=None):
    model = HiddenMarkovModel(name=name)
    s = {}

    # add states
    i0 = State(DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}), name=f"{name}:I0")

    model.add_state(i0)

    s[i0.name] = i0

    for c in range(len(target)):
        dc = State(None, name=f"{name}:D{c + 1}")

        mc = State(DiscreteDistribution({
            "A": 0.94 if target[c] == 'A' else 0.02,
            "C": 0.94 if target[c] == 'C' else 0.02,
            "G": 0.94 if target[c] == 'G' else 0.02,
            "T": 0.94 if target[c] == 'T' else 0.02,
        }), name=f"{name}:M{c + 1}")

        ic = State(DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}), name=f"{name}:I{c + 1}")

        model.add_states([mc, ic, dc])

        s[dc.name] = dc
        s[mc.name] = mc
        s[ic.name] = ic

    # add transitions
    model.add_transition(model.start, s[f'{name}:I0'], 0.05)
    model.add_transition(model.start, s[f'{name}:D1'], 0.05)
    model.add_transition(model.start, s[f'{name}:M1'], 0.90)

    model.add_transition(s[f'{name}:I0'], s[f'{name}:I0'], 0.70)
    model.add_transition(s[f'{name}:I0'], s[f'{name}:D1'], 0.15)
    model.add_transition(s[f'{name}:I0'], s[f'{name}:M1'], 0.15)

    for c in range(1, len(target)):
        model.add_transition(s[f'{name}:D{c}'], s[f'{name}:D{c + 1}'], 0.15)
        model.add_transition(s[f'{name}:D{c}'], s[f'{name}:I{c}'], 0.70)
        model.add_transition(s[f'{name}:D{c}'], s[f'{name}:M{c + 1}'], 0.15)

        model.add_transition(s[f'{name}:I{c}'], s[f'{name}:D{c + 1}'], 0.15)
        model.add_transition(s[f'{name}:I{c}'], s[f'{name}:I{c}'], 0.15)
        model.add_transition(s[f'{name}:I{c}'], s[f'{name}:M{c + 1}'], 0.70)

        model.add_transition(s[f'{name}:M{c}'], s[f'{name}:D{c + 1}'], 0.05)
        model.add_transition(s[f'{name}:M{c}'], s[f'{name}:I{c}'], 0.05)
        model.add_transition(s[f'{name}:M{c}'], s[f'{name}:M{c + 1}'], 0.90)

    model.add_transition(s[f'{name}:D{len(target)}'], s[f'{name}:I{len(target)}'], 0.70)
    model.add_transition(s[f'{name}:D{len(target)}'], model.end, 0.30)

    model.add_transition(s[f'{name}:I{len(target)}'], s[f'{name}:I{len(target)}'], 0.15)
    model.add_transition(s[f'{name}:I{len(target)}'], model.end, 0.85)

    model.add_transition(s[f'{name}:M{len(target)}'], s[f'{name}:I{len(target)}'], 0.90)
    model.add_transition(s[f'{name}:M{len(target)}'], model.end, 0.10)

    model.bake(merge='None')

    return model


def make_random_repeat_model(name='random'):
    model = HiddenMarkovModel(name=name)

    # add states
    ri = State(DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}), name=f"{name}:RI")
    rda = State(None, name=f"{name}:RDA")
    rdb = State(None, name=f"{name}:RDB")

    model.add_states([ri, rda, rdb])

    # add transitions
    model.add_transition(model.start, rda, 0.5)
    model.add_transition(model.start, ri, 0.5)

    model.add_transition(ri, ri, 0.8)
    model.add_transition(ri, rda, 0.10)
    model.add_transition(ri, model.end, 0.10)

    model.add_transition(rdb, ri, 0.5)
    model.add_transition(rdb, model.end, 0.5)

    model.bake(merge='None')

    return model


def build_default_model():
    adapters = {
        "10x_Adapter": "TCTACACGACGCTCTTCCGATCT",
        "5p_TSO": "TTTCTTATATGGG",
        "Poly_A": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
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
        "P": "ATGAAGTGGCTCGAGA"
    }

    direct_connections = {
        "Poly_A": ["3p_Adapter"],
        "3p_Adapter": ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"],
        "A": ["10x_Adapter"],
        "B": ["10x_Adapter"],
        "C": ["10x_Adapter"],
        "D": ["10x_Adapter"],
        "E": ["10x_Adapter"],
        "F": ["10x_Adapter"],
        "G": ["10x_Adapter"],
        "H": ["10x_Adapter"],
        "I": ["10x_Adapter"],
        "J": ["10x_Adapter"],
        "K": ["10x_Adapter"],
        "L": ["10x_Adapter"],
        "M": ["10x_Adapter"],
        "N": ["10x_Adapter"],
        "O": ["10x_Adapter"],
        "P": ["10x_Adapter"]
    }

    full_model = make_random_repeat_model()
    for k in adapters:
        full_model.add_model(make_global_alignment_model(adapters[k], k))

    full_model.bake(merge='None')

    # dictionary of model starting states, random start, and random end
    starts = {}
    rda = None
    rdb = None
    for s in full_model.states:
        if "-start" in s.name and "random" not in s.name:
            starts[re.sub("-start", "", s.name)] = s
        elif "random:RDA" in s.name:
            rda = s
        elif "random:RDB" in s.name:
            rdb = s

    # link rda to starts
    for sname in starts:
        full_model.add_transition(rda, starts[sname], 1.0 / len(starts))

    # link up ending states according to our direct connections dictionary
    for s in full_model.states:
        m = re.match("^(\w+):([MID])(\d+)", s.name)
        if m is not None and int(m.group(3)) == len(adapters[m.group(1)]):
            sname = m.group(1)

            if sname in direct_connections:
                p = 1.0 / (10 * len(direct_connections[sname]))
                full_model.add_transition(s, rdb, p)

                for dcname in direct_connections[sname]:
                    full_model.add_transition(s, starts[dcname], (1.0 - p) / len(direct_connections[sname]))
            else:
                full_model.add_transition(s, rdb, 0.5)

    full_model.bake()

    return full_model


def annotate(full_model, seq):
    logp, path = full_model.viterbi(seq)

    ppath = []
    for p, (idx, state) in enumerate(path[1:-1]):
        if "start" not in state.name and ":RD" not in state.name and ":D" not in state.name:
            ppath.append(f'{re.split(":", state.name)[0]}')

    return logp, ppath


# IUPAC RC's from: http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
# and https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
RC_BASE_MAP = {"N": "N", "A": "T", "T": "A", "G": "C", "C": "G", "Y": "R", "R": "Y", "S": "S", "W": "W", "K": "M",
               "M": "K", "B": "V", "V": "B", "D": "H", "H": "D", "n": "n", "a": "t", "t": "a", "g": "c", "c": "g",
               "y": "r", "r": "y", "s": "s", "w": "w", "k": "m", "m": "k", "b": "v", "v": "b", "d": "h", "h": "d"}


def reverse_complement(base_string):
    """
    Reverse complements the given base_string.
    :param base_string: String of bases to be reverse-complemented.
    :return: The reverse complement of the given base string.
    """

    return ''.join(map(lambda b: RC_BASE_MAP[b], base_string[::-1]))

