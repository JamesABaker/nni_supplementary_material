from __future__ import division
import scipy
from scipy import stats
import numpy as np
import pylab as P
import random

print "Positive residue statistics.\n"
print "File, Single-pass (SP) helices, SP Positive residues, SP t-statistic, SP P-value, Multi-pass (MP) ids, MP helices total, MP average helices per id, MP std of average helices per id, MP Positive residues, MP t-statstic, MP P-value"

list_of_files = [
    "TopDB_10_flanklength.csv",
    "UniHuman_10_flanklength.csv",
    "UniER_10_flanklength.csv",
    "UniGolgi_10_flanklength.csv",
    "UniPM_10_flanklength.csv",
    "UniCress_10_flanklength.csv",
    "UniFungi_10_flanklength.csv",
    "UniBacilli_10_flanklength.csv",
    "UniEcoli_10_flanklength.csv",
    "UniArch_10_flanklength.csv",
    ]

for file in list_of_files:

    list_of_tmh_segments = []
    list_of_singlepass_tmh_segments = []
    list_of_multipass_tmh_segments = []
    list_of_total_residues_multipass_inside_positive = []
    list_of_total_residues_multipass_inside_negative = []
    list_of_total_residues_multipass_outside_positive = []
    list_of_total_residues_multipass_outside_negative = []

    list_of_total_residues_singlepass_outside_positive = []
    list_of_total_residues_singlepass_inside_positive = []
    list_of_total_residues_singlepass_inside_negative = []
    list_of_total_residues_singlepass_outside_negative = []
    results = []
    list_of_multipass_ids = []
    list_of_multipass_helix_counts = []
    single_total_inside_flanks = ""
    single_total_outside_flanks = ""
    multi_total_inside_flanks = ""
    multi_total_outside_flanks = ""

    with open(file) as inputfile:
        for line in inputfile:
            results.append(line.strip().split(','))

    maximum_tmd_length = 0
    total_positive_single_flanks = 0
    total_negative_single_flanks = 0
    total_positive_multi_flanks = 0
    total_negative_multi_flanks = 0

    for entry in results:
        if entry == results[0]:
            pass
        else:
            name = entry[0]
            id = entry[1]
            n_terminal_start = entry[2]
            tmh_start_location = entry[3]
            tmh_end_location = entry[4]
            sequence = entry[5]
            tmh_sequence = entry[6]
            n_flank_sequence = entry[7]
            c_flank_sequence = entry[8]
            tmh_number = entry[9]
            total_tmd_count = entry[10]
            correction_number = 0

            if "Outside" in str(n_terminal_start):
                tmh_unaltered_sequence = str(tmh_sequence)
                tmh_reversed_sequence = tmh_unaltered_sequence[::-1]
                tmh_segment = tmh_reversed_sequence
                c_unaltered_sequence = str(c_flank_sequence)
                c_reversed_sequence = c_unaltered_sequence[::-1]
                inside_segment = c_reversed_sequence
                n_unaltered_sequence = str(n_flank_sequence)
                n_reversed_sequence = n_unaltered_sequence[::-1]
                outside_segment = n_reversed_sequence

            if "Inside" in str(n_terminal_start):
                tmh_unaltered_sequence = str(n_flank_sequence)
                tmh_segment = tmh_unaltered_sequence
                c_unaltered_sequence = str(c_flank_sequence)
                outside_segment = c_unaltered_sequence
                n_unaltered_sequence = str(n_flank_sequence)
                inside_segment = n_unaltered_sequence

            total_negative_residues_inside = inside_segment.count(
                "D") + inside_segment.count("E")
            total_negative_residues_outside = outside_segment.count(
                "D") + outside_segment.count("E")
            total_positive_residues_inside = inside_segment.count(
                "K") + inside_segment.count("R")
            total_positive_residues_outside = outside_segment.count(
                "K") + outside_segment.count("R")

            if int(total_tmd_count) == 1:
                list_of_singlepass_tmh_segments.append(tmh_segment)
                list_of_total_residues_singlepass_outside_positive.append(
                    total_positive_residues_outside)
                list_of_total_residues_singlepass_inside_positive.append(
                    total_positive_residues_inside)
                list_of_total_residues_singlepass_inside_negative.append(
                    total_negative_residues_inside)
                list_of_total_residues_singlepass_outside_negative.append(
                    total_negative_residues_outside)

                total_positive_single_flanks = total_positive_residues_inside + \
                    total_positive_residues_outside + total_positive_single_flanks
                total_negative_single_flanks = total_negative_residues_inside + \
                    total_negative_residues_outside + total_negative_single_flanks

                single_total_inside_flanks = single_total_inside_flanks + inside_segment
                single_total_outside_flanks = single_total_outside_flanks + outside_segment

            elif int(total_tmd_count) > 1:
                list_of_multipass_tmh_segments.append(tmh_segment)
                list_of_total_residues_multipass_outside_positive.append(
                    total_positive_residues_outside)
                list_of_total_residues_multipass_inside_positive.append(
                    total_positive_residues_inside)
                list_of_total_residues_multipass_inside_negative.append(
                    total_negative_residues_inside)
                list_of_total_residues_multipass_outside_negative.append(
                    total_negative_residues_outside)

                total_positive_multi_flanks = total_positive_residues_inside + \
                    total_positive_residues_outside + total_positive_multi_flanks
                total_negative_multi_flanks = total_negative_residues_inside + \
                    total_negative_residues_outside + total_negative_multi_flanks

                multi_total_inside_flanks = multi_total_inside_flanks + inside_segment
                multi_total_outside_flanks = multi_total_outside_flanks + outside_segment

                if id not in list_of_multipass_ids:
                    list_of_multipass_ids.append(id)
                    list_of_multipass_helix_counts.append(int(total_tmd_count))

    ttest_singlepass_positive = scipy.stats.ttest_ind(
        list_of_total_residues_singlepass_inside_positive, list_of_total_residues_singlepass_outside_positive)
    ttest_singlepass_negative = scipy.stats.ttest_ind(
        list_of_total_residues_singlepass_inside_negative, list_of_total_residues_singlepass_outside_negative)

    ttest_multipass_positive = scipy.stats.ttest_ind(
        list_of_total_residues_multipass_inside_positive, list_of_total_residues_multipass_outside_positive)
    ttest_multipass_negative = scipy.stats.ttest_ind(
        list_of_total_residues_multipass_inside_negative, list_of_total_residues_multipass_outside_negative)

    ks_singlepass_positive_inside = scipy.stats.kstest(
        list_of_total_residues_singlepass_inside_positive, 'norm', alternative='greater')
    ks_singlepass_positive_outside = scipy.stats.kstest(
        list_of_total_residues_singlepass_outside_positive, 'norm', alternative='greater')
    ks_singlepass_negative_inside = scipy.stats.kstest(
        list_of_total_residues_singlepass_inside_negative, 'norm', alternative='greater')
    ks_singlepass_negative_outside = scipy.stats.kstest(
        list_of_total_residues_singlepass_outside_negative, 'norm', alternative='greater')

    ks_multipass_positive_inside = scipy.stats.kstest(
        list_of_total_residues_multipass_inside_positive, 'norm', alternative='greater')
    ks_multipass_positive_outside = scipy.stats.kstest(
        list_of_total_residues_multipass_outside_positive, 'norm', alternative='greater')
    ks_multipass_negative_inside = scipy.stats.kstest(
        list_of_total_residues_multipass_inside_negative, 'norm', alternative='greater')
    ks_multipass_negative_outside = scipy.stats.kstest(
        list_of_total_residues_multipass_outside_negative, 'norm', alternative='greater')


    average_helix_count = np.mean(list_of_multipass_helix_counts)
    std_of_helix_count = np.std(list_of_multipass_helix_counts)


    print file,",", len(list_of_singlepass_tmh_segments),",", total_positive_single_flanks,",", ttest_singlepass_positive[0],",", ttest_singlepass_positive[1],",", len(list_of_multipass_ids),",", len(list_of_multipass_tmh_segments),",", average_helix_count,",", std_of_helix_count,",", total_positive_multi_flanks,",", ttest_multipass_positive[0],",", ttest_multipass_positive[1]
