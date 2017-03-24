from __future__ import division
import scipy
from scipy import stats
import numpy as np
import pylab as P
import random

print "Negative residue statistics.\n"
print "File, Single-pass (singlepass) helices, singlepass negative residues inside, singlepass negative residues outside, singlepass t-statistic, singlepass P-value, Multi-pass (multipass) ids, multipass helices total, multipass average helices per id, multipass std of average helices per id, multipass negative residues inside, multipass negative residues outside,, multipass t-statstic, multipass P-value"


list_of_files = [
    "TopDB_20_flanklength_flankclashTrue.csv",
    "UniHuman_20_flanklength_flankclashTrue.csv",
    "UniER_20_flanklength_flankclashTrue.csv",
    "UniGolgi_20_flanklength_flankclashTrue.csv",
    "UniPM_20_flanklength_flankclashTrue.csv",
    "UniCress_20_flanklength_flankclashTrue.csv",
    "UniFungi_20_flanklength_flankclashTrue.csv",
    "UniBacilli_20_flanklength_flankclashTrue.csv",
    "UniEcoli_20_flanklength_flankclashTrue.csv",
    "UniArch_20_flanklength_flankclashTrue.csv",
]


for file in list_of_files:

    list_of_tmh_sequence_to_uses = []
    list_of_singlepass_tmh_sequence_to_uses = []
    list_of_multipass_tmh_sequence_to_uses = []
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

    total_e_residues_inside_singlepass = 0
    total_e_residues_outside_singlepass = 0
    total_d_residues_inside_singlepass = 0
    total_d_residues_outside_singlepass = 0

    total_e_residues_inside_multipass = 0
    total_e_residues_outside_multipass = 0
    total_d_residues_inside_multipass = 0
    total_d_residues_outside_multipass = 0

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

            tmh_sequence_to_use = n_flank_sequence + tmh_sequence + c_flank_sequence

            def inside_limit_positions(tmh_string, inside_flank_sequence_string):
                '''
                This returns the positions to obtain the inside flank sensitively to avoid scrambling the slice with a below 0 slice index.
                '''
                # The half way point can be quickly found by adding length of
                # the "left-hand"/inside flank to half the length of the TMH.
                # The +1 insures that when dividing by 2 that odd and even
                # total length numbers both end up at "position 0" rather than
                # -1 or 0.
                real_lower_value = int(
                    ((len(tmh_string) + 1) / 2) + len(inside_flank_sequence_string) - 20)
                if real_lower_value < 0:
                    lower_value = 0
                elif real_lower_value >= 0:
                    lower_value = real_lower_value

                real_upper_value = int(
                    ((len(tmh_string) + 1) / 2) + len(inside_flank_sequence_string) - 10)
                if real_upper_value < 0:
                    upper_value = 0
                elif real_upper_value >= 0:
                    upper_value = real_upper_value

                # Now these slices will be at LEAST [0:N]
                slice_values = [lower_value, upper_value]
                return slice_values

            if "Outside" in str(n_terminal_start):
                this_sequence = tmh_sequence_to_use[::-1]

                # The inside segment needs lower and upper bounds to be handled
                # sensitively since a negative value will cause
                # nothing/scrambled strings to be recorded. The outside end
                # doesn't really matter if the values are too high; the values
                # recorded will ascend to the maximum value and stop.
                values_of_inside_limit_positions = inside_limit_positions(
                    tmh_sequence, c_flank_sequence)
                # print
                # values_of_inside_limit_positions[0],values_of_inside_limit_positions[1]
                inside_segment = this_sequence[values_of_inside_limit_positions[
                    0]:values_of_inside_limit_positions[1]]

                # print int(((len(tmh_sequence)+1) / 2) + len(c_flank_sequence)
                # + 10), int(((len(tmh_sequence)+1) / 2) +
                # len(c_flank_sequence) + 20)
                outside_segment = this_sequence[int(((len(tmh_sequence) + 1) / 2) + len(
                    c_flank_sequence) + 10):int(((len(tmh_sequence) + 1) / 2) + len(c_flank_sequence) + 20)]

            if "Inside" in str(n_terminal_start):
                this_sequence = tmh_sequence_to_use

                values_of_inside_limit_positions = inside_limit_positions(
                    tmh_sequence, n_flank_sequence)
                inside_segment = this_sequence[values_of_inside_limit_positions[
                    0]:values_of_inside_limit_positions[1]]

                outside_segment = this_sequence[int(((len(tmh_sequence) + 1) / 2) + len(
                    n_flank_sequence) + 10):int(((len(tmh_sequence) + 1) / 2) + len(n_flank_sequence) + 20)]

            total_negative_residues_inside = inside_segment.count(
                "D") + inside_segment.count("E")
            total_negative_residues_outside = outside_segment.count(
                "D") + outside_segment.count("E")
            total_positive_residues_inside = inside_segment.count(
                "K") + inside_segment.count("R")
            total_positive_residues_outside = outside_segment.count(
                "K") + outside_segment.count("R")

            if int(total_tmd_count) == 1:
                list_of_singlepass_tmh_sequence_to_uses.append(
                    tmh_sequence_to_use)
                list_of_total_residues_singlepass_outside_positive.append(
                    total_positive_residues_outside)
                list_of_total_residues_singlepass_inside_positive.append(
                    total_positive_residues_inside)
                list_of_total_residues_singlepass_inside_negative.append(
                    total_negative_residues_inside)
                list_of_total_residues_singlepass_outside_negative.append(
                    total_negative_residues_outside)

            elif int(total_tmd_count) > 1:
                list_of_multipass_tmh_sequence_to_uses.append(
                    tmh_sequence_to_use)
                list_of_total_residues_multipass_outside_positive.append(
                    total_positive_residues_outside)
                list_of_total_residues_multipass_inside_positive.append(
                    total_positive_residues_inside)
                list_of_total_residues_multipass_inside_negative.append(
                    total_negative_residues_inside)
                list_of_total_residues_multipass_outside_negative.append(
                    total_negative_residues_outside)

    stat_test_singlepass_positive = scipy.stats.kruskal(
        list_of_total_residues_singlepass_inside_positive, list_of_total_residues_singlepass_outside_positive)
    stat_test_singlepass_negative = scipy.stats.kruskal(
        list_of_total_residues_singlepass_inside_negative, list_of_total_residues_singlepass_outside_negative)

    stat_test_multipass_positive = scipy.stats.kruskal(
        list_of_total_residues_multipass_inside_positive, list_of_total_residues_multipass_outside_positive)
    stat_test_multipass_negative = scipy.stats.kruskal(
        list_of_total_residues_multipass_inside_negative, list_of_total_residues_multipass_outside_negative)

    singlepass_negative_residues_inside = sum(
        list_of_total_residues_singlepass_inside_negative)
    singlepass_negative_residues_outside = sum(
        list_of_total_residues_singlepass_outside_negative)

    multipass_negative_residues_inside = sum(
        list_of_total_residues_multipass_inside_negative)
    multipass_negative_residues_outside = sum(
        list_of_total_residues_multipass_outside_negative)

    print file, ",", len(list_of_singlepass_tmh_sequence_to_uses), ",", singlepass_negative_residues_inside, ",", singlepass_negative_residues_outside, ",", stat_test_singlepass_negative[0], ",", stat_test_singlepass_negative[1], ",", multipass_negative_residues_inside, ",", multipass_negative_residues_outside, ",", stat_test_multipass_negative[0], ",", stat_test_multipass_negative[1]
