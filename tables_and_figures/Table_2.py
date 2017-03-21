from __future__ import division
import scipy
from scipy import stats
import numpy as np
import pylab as P
import random

print "Negative residue statistics.\n"
print "File, Single-pass (SP) helices, SP negative residues inside, SP negative residues outside, SP t-statistic, SP P-value, Multi-pass (MP) IDs, MP helices total, MP average helices per ID, MP std of average helices per ID, MP negative residues inside, MP negative residues outside,, MP t-statstic, MP P-value"


list_of_files = [
    "TopDB_10_flanklength_flankclashFalse.csv",
    "UniHuman_10_flanklength_flankclashFalse.csv",
    "UniER_10_flanklength_flankclashFalse.csv",
    "UniGolgi_10_flanklength_flankclashFalse.csv",
    "UniPM_10_flanklength_flankclashFalse.csv",
    "UniCress_10_flanklength_flankclashFalse.csv",
    "UniFungi_10_flanklength_flankclashFalse.csv",
    "UniBacilli_10_flanklength_flankclashFalse.csv",
    "UniEcoli_10_flanklength_flankclashFalse.csv",
    "UniArch_10_flanklength_flankclashFalse.csv",

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



    "top_all_10_flanklength_flankclashTrue_only_half_flanks.csv",
    "UniHuman_10_flanklength_flankclashTrue_only_half_flanks.csv",


    "UniER_10_flanklength_flankclashFalse_only_half_flanks.csv",
    "UniER_10_flanklength_flankclashTrue_only_half_flanks.csv",
    "UniGolgi_10_flanklength_flankclashFalse_only_half_flanks.csv",
    "UniGolgi_10_flanklength_flankclashTrue_only_half_flanks.csv",
    "UniPM_10_flanklength_flankclashFalse_only_half_flanks.csv",
    "UniPM_10_flanklength_flankclashTrue_only_half_flanks.csv",
    "UniCress_10_flanklength_flankclashFalse_only_half_flanks.csv",
    "UniCress_10_flanklength_flankclashTrue_only_half_flanks.csv",
    "UniFungi_10_flanklength_flankclashFalse_only_half_flanks.csv",
    "UniFungi_10_flanklength_flankclashTrue_only_half_flanks.csv",
    "UniBacilli_10_flanklength_flankclashFalse_only_half_flanks.csv",
    "UniBacilli_10_flanklength_flankclashTrue_only_half_flanks.csv",
    "UniEcoli_10_flanklength_flankclashFalse_only_half_flanks.csv",
    "UniEcoli_10_flanklength_flankclashTrue_only_half_flanks.csv",
    "UniArch_10_flanklength_flankclashFalse_only_half_flanks.csv",
    "UniArch_10_flanklength_flankclashTrue_only_half_flanks.csv",
    "UniHuman_10_flanklength_flankclashFalse_only_full_flanks.csv",
    "UniHuman_10_flanklength_flankclashTrue_only_full_flanks.csv",
    "UniER_10_flanklength_flankclashFalse_only_full_flanks.csv",
    "UniER_10_flanklength_flankclashTrue_only_full_flanks.csv",
    "UniGolgi_10_flanklength_flankclashFalse_only_full_flanks.csv",
    "UniGolgi_10_flanklength_flankclashTrue_only_full_flanks.csv",
    "UniPM_10_flanklength_flankclashFalse_only_full_flanks.csv",
    "UniPM_10_flanklength_flankclashTrue_only_full_flanks.csv",
    "UniCress_10_flanklength_flankclashFalse_only_full_flanks.csv",
    "UniCress_10_flanklength_flankclashTrue_only_full_flanks.csv",
    "UniFungi_10_flanklength_flankclashFalse_only_full_flanks.csv",
    "UniFungi_10_flanklength_flankclashTrue_only_full_flanks.csv",
    "UniBacilli_10_flanklength_flankclashFalse_only_full_flanks.csv",
    "UniBacilli_10_flanklength_flankclashTrue_only_full_flanks.csv",
    "UniEcoli_10_flanklength_flankclashFalse_only_full_flanks.csv",
    "UniEcoli_10_flanklength_flankclashTrue_only_full_flanks.csv",
    "UniArch_10_flanklength_flankclashFalse_only_full_flanks.csv",
    "UniArch_10_flanklength_flankclashTrue_only_full_flanks.csv",
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
        temp = inputfile.read().splitlines()
        for line in temp:
            results.append(line.strip().split(','))
    # print len(results)

    maximum_tmd_length = 0
    total_positive_single_flanks = 0
    total_negative_single_flanks = 0
    total_positive_multi_flanks = 0
    total_negative_multi_flanks = 0
    for entry in results:
        # print entry
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
            tmh_number = int(entry[9])
            total_tmd_count = int(entry[10])
            correction_number = 0

            if "Outside" in str(n_terminal_start):
                tmh_unaltered_sequence = str(tmh_sequence)
                tmh_reversed_sequence = tmh_unaltered_sequence[::-1]
                tmh_segment = tmh_reversed_sequence
                C_unaltered_sequence = str(C_flank_sequence)
                C_reversed_sequence = C_unaltered_sequence[::-1]
                inside_segment = C_reversed_sequence
                N_unaltered_sequence = str(N_flank_sequence)
                N_reversed_sequence = N_unaltered_sequence[::-1]
                outside_segment = N_reversed_sequence

            if "Inside" in str(n_terminal_start):
                tmh_unaltered_sequence = str(N_flank_sequence)
                tmh_segment = tmh_unaltered_sequence
                C_unaltered_sequence = str(C_flank_sequence)
                outside_segment = C_unaltered_sequence
                N_unaltered_sequence = str(N_flank_sequence)
                inside_segment = N_unaltered_sequence

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

    if len(list_of_total_residues_singlepass_inside_negative) > 0 and len(list_of_total_residues_singlepass_outside_negative) > 0:
        ttest_singlepass_positive = scipy.stats.kruskal(
            list_of_total_residues_singlepass_inside_positive, list_of_total_residues_singlepass_outside_positive)
        ttest_singlepass_negative = scipy.stats.kruskal(
            list_of_total_residues_singlepass_inside_negative, list_of_total_residues_singlepass_outside_negative)

        ttest_multipass_positive = scipy.stats.kruskal(
            list_of_total_residues_multipass_inside_positive, list_of_total_residues_multipass_outside_positive)
        ttest_multipass_negative = scipy.stats.kruskal(
            list_of_total_residues_multipass_inside_negative, list_of_total_residues_multipass_outside_negative)

        sp_negative_residues_inside = sum(
            list_of_total_residues_singlepass_inside_negative)
        sp_negative_residues_outside = sum(
            list_of_total_residues_singlepass_outside_negative)

        mp_negative_residues_inside = sum(
            list_of_total_residues_multipass_inside_negative)
        mp_negative_residues_outside = sum(
            list_of_total_residues_multipass_outside_negative)

        print file, ",", len(list_of_singlepass_tmh_segments), ",", sp_negative_residues_inside, ",", sp_negative_residues_outside, ",", ttest_singlepass_negative[0], ",", ttest_singlepass_negative[1], ",", len(list_of_multipass_tmh_segments), ",", mp_negative_residues_inside, ",", mp_negative_residues_outside, ",", ttest_multipass_negative[0], ",", ttest_multipass_negative[1]
    else:
        print file, "HAD NO VIABLE FLANK PAIRS FOR COMPARISON"
