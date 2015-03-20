#!/usr/bin/env python
# -*- coding: utf8 -*-

from __future__ import division

from math import log, exp
from decimal import Decimal
import sys

from Bio import SeqIO
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


def pretty_scientific_p(p):
    return '{0:.4e}'.format(p)


def pretty_e_p(p):
    return pretty_p(exp(p))


def pretty_scientific_e_p(p):
    return pretty_scientific_p(exp(p))


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
            transition_statistics = [[log(value) - log(state_statistics[k]) for value in row] for k, row in
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
                print('\t'.join(pretty_e_p(value) for value in row))

        report()
    print(gc_rich_segments)


CONVERGENCE_LIKELYHOOD_INCREASE = Decimal(0.1)


def run_baum_welch(emission_p, sequence, start_p, transition_p):
    rounds = 0
    length = len(sequence)
    log_likelyhood = -10050000000000
    while True:
        rounds += 1
        alpha = [[0 for _ in STATES] for _ in xrange(length)]
        beta = [[0 for _ in STATES] for _ in xrange(length)]

        for state in STATES:
            alpha[0][state] = start_p[state] + emission_p[state][sequence[0]]

        for i in xrange(1, len(sequence)):
            nucleotide = sequence[i]
            for j in STATES:
                alpha[i][j] = emission_p[j][nucleotide] + iter_slog(
                    alpha[i - 1][state] + transition_p[state][j] for state in STATES)

        for state in STATES:
            beta[length - 1][state] = 1

        for i in (length - j - 1 for j in xrange(1, len(sequence))):
            nucleotide = sequence[i + 1]
            for j in STATES:
                beta[i][j] = iter_slog(
                    beta[i + 1][state] + transition_p[j][state] + emission_p[j][nucleotide] for state in STATES)

        gamma = [[0 for _ in STATES] for _ in xrange(length)]
        ksi = [[[0 for _ in STATES] for _ in STATES] for _ in xrange(length - 1)]
        sums = [0 for _ in xrange(length)]
        for i in xrange(length):
            for j in STATES:
                gamma[i][j] = alpha[i][j] + beta[i][j]
            sums[i] = iter_slog(gamma[i])
            for j in STATES:
                gamma[i][j] -= sums[i]

        for i in xrange(length - 1):
            nucleotide = sequence[i + 1]
            for prev in STATES:
                for next in STATES:
                    ksi[i][prev][next] = alpha[i][prev] + beta[i + 1][next] + transition_p[prev][next] + \
                                         emission_p[next][nucleotide] - sums[i]

        start_p = gamma[0]
        for i in STATES:
            for j in STATES:
                transition_p[i][j] = iter_slog(ksi[k][i][j] for k in xrange(length - 1)) - iter_slog(
                    gamma[k][i] for k in xrange(length - 1))

        for state in STATES:
            for nucleotide in ALPHABET:
                emission_p[state][nucleotide] = iter_slog(
                    gamma[k][state] for k in xrange(length) if sequence[k] == nucleotide) - iter_slog(
                    gamma[k][state] for k in xrange(length))

        new_log_likelyhood = iter_slog(alpha[length - 1])
        if new_log_likelyhood < log_likelyhood:
            print('BAD STUFF')
            break
        else:
            print(new_log_likelyhood)
            if new_log_likelyhood - log_likelyhood < CONVERGENCE_LIKELYHOOD_INCREASE and rounds > 1:
                print('CONVERGED!')
                print('Took us {0} rounds'.format(rounds))
                print('Transition probabilities:')
                for row in transition_p:
                    print('\t'.join(pretty_scientific_e_p(value) for value in row))
                print('Emission probabilities:')
                for row in emission_p:
                    print('\t'.join('{0}: {1}'.format(key, pretty_scientific_e_p(value))
                                    for key, value in row.iteritems()))
                break
            log_likelyhood = new_log_likelyhood


def main():
    if len(sys.argv) == 1 or len(sys.argv) > 3:
        print('Usage:\n\tlab2.py <input_fasta> [<output>]')
        exit(1)

    filename = sys.argv[1]
    if len(sys.argv) == 3:
        sys.stdout = open(sys.argv[2], 'w')
    start_p = list(map(log, [0.996, 0.004]))
    transition_p = [list(map(log, [0.999, 0.001])),
                    list(map(log, [0.01, 0.99]))]
    emission_p = [{'a': log(0.291), 'c': log(0.209), 'g': log(0.209), 't': log(0.291)},
                  {'a': log(0.169), 'c': log(0.331), 'g': log(0.331), 't': log(0.169)}]

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