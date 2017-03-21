from __future__ import division
import scipy
from scipy import stats
import numpy as np
import pylab as P
import random

print "Negative residue statistics.\n"
print "File, Single-pass (SP) helices, SP negative residues inside, SP negative residues outside, SP t-statistic, SP P-value, Multi-pass (MP) IDs, MP helices total, MP average helices per ID, MP std of average helices per ID, MP negative residues inside, MP negative residues outside,, MP t-statstic, MP P-value"


list_of_files = [
    "TopDB_10_flanklength_flankclashTrue.csv",
    "UniHuman_10_flanklength_flankclashTrue.csv",
    "UniER_10_flanklength_flankclashTrue.csv",
    "UniGolgi_10_flanklength_flankclashTrue.csv",
    "UniPM_10_flanklength_flankclashTrue.csv",
    "UniCress_10_flanklength_flankclashTrue.csv",
    "UniFungi_10_flanklength_flankclashTrue.csv",
    "UniBacilli_10_flanklength_flankclashTrue.csv",
    "UniEcoli_10_flanklength_flankclashTrue.csv",
    "UniArch_10_flanklength_flankclashTrue.csv",
]


for file in list_of_files:

    list_of_tmh_sequence_to_uses = []
    list_of_singlepass_tmh_sequence_to_use = []
    list_of_multipass_tmh_sequence_to_use = []
    list_of_total_residues_multipass_inside_positive = []
    list_of_total_residues_multipass_inside_negative = []
    list_of_total_residues_multipass_outside_positive = []
    list_of_total_residues_multipass_outside_negative = []
    list_of_multipass_tmh_negative_residues = []

    list_of_total_residues_singlepass_outside_positive = []
    list_of_total_residues_singlepass_inside_positive = []
    list_of_total_residues_singlepass_inside_negative = []
    list_of_total_residues_singlepass_outside_negative = []
    list_of_singlepass_tmh_negative_residues = []

    results = []
    list_of_multipass_IDs = []
    list_of_multipass_helix_counts = []
    single_total_inside_flanks = ""
    single_total_outside_flanks = ""
    multi_total_inside_flanks = ""
    multi_total_outside_flanks = ""

    total_E_residues_inside_singlepass = 0
    total_E_residues_outside_singlepass = 0
    total_D_residues_inside_singlepass = 0
    total_D_residues_outside_singlepass = 0

    total_E_residues_inside_multipass = 0
    total_E_residues_outside_multipass = 0
    total_D_residues_inside_multipass = 0
    total_D_residues_outside_multipass = 0

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
            N_flank_sequence = entry[7]
            C_flank_sequence = entry[8]
            tmh_number = entry[9]
            total_tmd_count = entry[10]
            correction_number = 0

            tmh_sequence_to_use = N_flank_sequence + tmh_sequence + C_flank_sequence

            if "Outside" in str(n_terminal_start):
                alignment_compensation_value = len(
                    C_flank_sequence) - len(N_flank_sequence)
                this_sequence = tmh_sequence_to_use[::-1]
                inside_segment = this_sequence[int((len(this_sequence) / 2) - 20) + alignment_compensation_value:int(
                    (len(this_sequence) / 2) - 10) + alignment_compensation_value]
                outside_segment = this_sequence[int((len(this_sequence) / 2) + 10) + alignment_compensation_value:int(
                    (len(this_sequence) / 2) + 20) + alignment_compensation_value]
                tmh_core = this_sequence[int((len(this_sequence) / 2) - 10) + alignment_compensation_value:int(
                    (len(this_sequence) / 2) + 10) + alignment_compensation_value]

            if "Inside" in str(n_terminal_start):
                alignment_compensation_value = len(
                    N_flank_sequence) - len(C_flank_sequence)
                this_sequence = tmh_sequence_to_use
                inside_segment = this_sequence[int((len(this_sequence) / 2) - 20) + alignment_compensation_value:int(
                    (len(this_sequence) / 2) - 10) + alignment_compensation_value]
                outside_segment = this_sequence[int((len(this_sequence) / 2) + 10) + alignment_compensation_value:int(
                    (len(this_sequence) / 2) + 20) + alignment_compensation_value]
                tmh_core = this_sequence[int((len(this_sequence) / 2) - 10) + alignment_compensation_value:int(
                    (len(this_sequence) / 2) + 10) + alignment_compensation_value]

            total_negative_residues_inside = inside_segment.count(
                "D") + inside_segment.count("E")
            total_negative_residues_outside = outside_segment.count(
                "D") + outside_segment.count("E")
            total_positive_residues_inside = inside_segment.count(
                "K") + inside_segment.count("R")
            total_positive_residues_outside = outside_segment.count(
                "K") + outside_segment.count("R")

            tmh_core_negative_residues = tmh_core.count(
                "D") + outside_segment.count("E")

            if int(total_tmd_count) == 1:
                list_of_singlepass_tmh_sequence_to_use.append(
                    tmh_sequence_to_use)

                list_of_total_residues_singlepass_inside_negative.append(
                    total_negative_residues_inside)
                list_of_total_residues_singlepass_outside_negative.append(
                    total_negative_residues_outside)

                list_of_singlepass_tmh_negative_residues.append(
                    tmh_core_negative_residues)

            elif int(total_tmd_count) > 1:
                list_of_multipass_tmh_sequence_to_use.append(
                    tmh_sequence_to_use)

                list_of_total_residues_multipass_inside_negative.append(
                    total_negative_residues_inside)
                list_of_total_residues_multipass_outside_negative.append(
                    total_negative_residues_outside)

                list_of_multipass_tmh_negative_residues.append(
                    tmh_core_negative_residues)

    stat_test_singlepass_positive = scipy.stats.kruskal(
        list_of_total_residues_singlepass_inside_positive, list_of_total_residues_singlepass_outside_positive)
    stat_test_singlepass_negative = scipy.stats.kruskal(
        list_of_total_residues_singlepass_inside_negative, list_of_total_residues_singlepass_outside_negative)

    stat_test_multipass_positive = scipy.stats.kruskal(
        list_of_total_residues_multipass_inside_positive, list_of_total_residues_multipass_outside_positive)
    stat_test_multipass_negative = scipy.stats.kruskal(
        list_of_total_residues_multipass_inside_negative, list_of_total_residues_multipass_outside_negative)

    sp_negative_residues_inside = sum(
        list_of_total_residues_singlepass_inside_negative)
    sp_negative_residues_outside = sum(
        list_of_total_residues_singlepass_outside_negative)

    mp_negative_residues_inside = sum(
        list_of_total_residues_multipass_inside_negative)
    mp_negative_residues_outside = sum(
        list_of_total_residues_multipass_outside_negative)
    # print "total negative residues in singlepass,",sum(list_of_singlepass_tmh_negative_residues)
    # print "total negative residues in
    # multipass,",sum(list_of_multipass_tmh_negative_residues)
    print file, ",", len(list_of_singlepass_tmh_sequence_to_use), ",", sp_negative_residues_inside, ",", sp_negative_residues_outside, ",", stat_test_singlepass_negative[0], ",", stat_test_singlepass_negative[1], ",", mp_negative_residues_inside, ",", mp_negative_residues_outside, ",", stat_test_multipass_negative[0], ",", stat_test_multipass_negative[1]
