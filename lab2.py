#!/usr/bin/env python
# -*- coding: utf8 -*-

from __future__ import division

from math import log, exp, log10
from Bio import SeqIO
from Bio.SeqUtils import GC
import sys
from blist import blist

ALPHABET = 'acgt'

paired = {'a': 't', 'c': 'g',
          'g': 'c', 't': 'a'}

AT_RICH = 0
GC_RICH = 1
STATES = [AT_RICH, GC_RICH]


def slog(logp, logq):
    # Richard Durbin
    # Biological sequence analysis: Probabilistic models of proteins and nucleic acids (1998), section 3.6
    if logq > logp:
        logp, logq = logq, logp
    return logp + log(1 + exp(logq - logp))


def iter_slog(l):
    it = iter(l)
    acc = next(it)
    for snd in it:
        acc = slog(acc, snd)
    return acc


def get_start_p(char, start, emission):
    row = []
    for state, state_emission in enumerate(emission):
        row.append(log(state_emission[char] * start[state]))
    return row


ROUNDS = 10


def pretty_p(p):
    return '{0:.4f}'.format(p)


def pretty_e_p(p):
    return pretty_p(10 ** p)


def run_viterbi(emission_p, sequence, start_p, transition_p):
    gc_rich_segments = []
    for r in xrange(ROUNDS):
        print('Round {0}'.format(r))
        viterbi = [[0 for _ in STATES]]
        state_sequence = [None for _ in STATES]
        new_state_sequence = [None for _ in STATES]

        for state in STATES:
            viterbi[0][state] = start_p[state] + emission_p[state][sequence[0]]
            state_sequence[state] = blist([state])

        for i in xrange(1, len(sequence)):
            nucleotide = sequence[i]
            viterbi.append([0 for _ in STATES])

            for state in STATES:
                (probability, prev_state) = max((viterbi[i - 1][prev_state] + transition_p[prev_state][state] +
                                                 emission_p[state][nucleotide], prev_state) for prev_state in STATES)
                viterbi[i][state] = probability
                new_list = blist(state_sequence[prev_state])
                new_list.append(state)
                new_state_sequence[state] = new_list
            state_sequence = new_state_sequence

        _, best_state = max((viterbi[len(viterbi) - 1][state], state) for state in STATES)
        # Report data and then change transitions

        def get_statistics(hidden_states_sequence):
            length = len(hidden_states_sequence)
            state = AT_RICH
            state_statistics = [0 for _ in STATES]
            segment_statistics = [0 for _ in STATES]
            segment_statistics[state] += 1
            transition_statistics = [[0 for _ in STATES] for _ in STATES]
            current_gc_rich_segment = [0, 0]
            gc_rich_segments = []
            for i, new_state in enumerate(hidden_states_sequence):
                state_statistics[new_state] += 1
                if state != new_state:
                    if new_state == GC_RICH:
                        current_gc_rich_segment[0] = i
                    else:
                        current_gc_rich_segment[1] = i
                        gc_rich_segments.append(tuple(current_gc_rich_segment))
                    transition_statistics[state][new_state] += 1
                    state = new_state
                    segment_statistics[new_state] += 1
                else:
                    transition_statistics[state][new_state] += 1
            transition_statistics = [[log10(value) - log10(state_statistics[k]) for value in row] for k, row in
                                     enumerate(transition_statistics)]
            return state_statistics, segment_statistics, transition_statistics, gc_rich_segments

        hidden_states_sequence = list(state_sequence[best_state])
        state_statistics, segment_statistics, transition_statistics, gc_rich_segments = \
            get_statistics(hidden_states_sequence)

        transition_p = transition_statistics

        def report():
            print('States: AT-rich = {0}\tGC-rich = {1}'.format(*state_statistics))
            print('Segments: AT-rich = {0}\tGC-rich = {1}'.format(*segment_statistics))
            for row in transition_statistics:
                print('{0}\t{1}'.format(*map(pretty_e_p, row)))
        report()
    print(gc_rich_segments)


def run_baum_welch(emission_p, sequence, start_p, transition_p):
    return


def main():
    if len(sys.argv) == 1 or len(sys.argv) > 3:
        print('Usage:\n\tlab2.py <input_fasta> [<output>]')
        exit(1)

    filename = sys.argv[1]
    if len(sys.argv) == 3:
        sys.stdout = open(sys.argv[2], 'w')
    start_p = list(map(log10, [0.996, 0.004]))
    transition_p = [list(map(log10, [0.999, 0.001])),
                    list(map(log10, [0.01, 0.99]))]
    emission_p = [{'a': log10(0.291), 'c': log10(0.209), 'g': log10(0.209), 't': log10(0.291)},
                  {'a': log10(0.169), 'c': log10(0.331), 'g': log10(0.331), 't': log10(0.169)}]

    with open(filename) as handle:
        sequence = SeqIO.read(handle, 'fasta').seq.lower()

    print('filename: {0}'.format(filename))
    with open(filename) as handle:
        print(handle.readline())

    print('====================')
    print('===== Viterbi ======')
    run_viterbi(emission_p, sequence, start_p, transition_p)

    print('====================')
    print('==== Baum-Welch ====')
    run_baum_welch(emission_p, sequence, start_p, transition_p)
    sys.stdout.close()


if __name__ == '__main__':
    main()