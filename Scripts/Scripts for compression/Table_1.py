from __future__ import division
import scipy
from scipy import stats
import numpy as np
import pylab as P
import random

list_of_files = [ #"TopDB_5_flanklength_flankclashFalse.csv",
                 #"TopDB_5_flanklength_flankclashTrue.csv",
                 "UniHuman_5_flanklength_flankclashFalse.csv",
                 "UniHuman_5_flanklength_flankclashTrue.csv",
                 "UniER_5_flanklength_flankclashFalse.csv",
                 "UniER_5_flanklength_flankclashTrue.csv",
                 "UniGolgi_5_flanklength_flankclashFalse.csv",
                 "UniGolgi_5_flanklength_flankclashTrue.csv",
                 "UniPM_5_flanklength_flankclashFalse.csv",
                 "UniPM_5_flanklength_flankclashTrue.csv",
                 "UniCress_5_flanklength_flankclashFalse.csv",
                 "UniCress_5_flanklength_flankclashTrue.csv",
                 "UniFungi_5_flanklength_flankclashFalse.csv",
                 "UniFungi_5_flanklength_flankclashTrue.csv",
                 "UniBacilli_5_flanklength_flankclashFalse.csv",
                 "UniBacilli_5_flanklength_flankclashTrue.csv",
                 "UniEcoli_5_flanklength_flankclashFalse.csv",
                 "UniEcoli_5_flanklength_flankclashTrue.csv",
                 "UniArch_5_flanklength_flankclashFalse.csv",
                 "UniArch_5_flanklength_flankclashTrue.csv",






                 ]

for file in list_of_files:

    single_pass_list_of_number_of_acidic_residues = []
    multi_pass_list_of_number_of_acidic_residues = []

    list_of_tmh_segments = []
    list_of_singlepass_tmh_segments = []
    list_of_multipass_tmh_segments = []
    results = []
    list_of_multipass_ids = []
    list_of_multipass_helix_counts = []

    single_total_leucines_leaflets = 0
    multi_total_leucines_leaflets = 0

    single_inside_leaflets = []
    single_outside_leaflets = []
    multi_inside_leaflets = []
    multi_outside_leaflets = []

    with open(file) as inputfile:
        for line in inputfile:
            results.append(line.strip().split(','))

    maximum_tmd_length = 0
    total_single_flanks = 0
    total_multi_flanks = 0

    for entry in results:
        if entry == results[0]:
            pass
        else:
            name = entry[0]
            id = entry[1]
            n_terminal_start = str(entry[2])
            tmh_start_location = int(entry[3])
            tmh_end_location = int(entry[4])
            sequence = str(entry[5])
            tmh_sequence = str(entry[6])
            n_flank_sequence = str(entry[7])
            c_flank_sequence = str(entry[8])
            tmh_number = int(entry[9])
            total_tmd_count = int(entry[10])
            correction_number = 0

            tmh_length = len(tmh_sequence)

            if int(total_tmd_count) == 1:
                acid_count_this_helix = tmh_sequence.count(
                    'E') + tmh_sequence.count('D')
                single_pass_list_of_number_of_acidic_residues.append(
                    acid_count_this_helix)

            elif int(total_tmd_count) > 1:
                acid_count_this_helix = tmh_sequence.count(
                    'E') + tmh_sequence.count('D')
                multi_pass_list_of_number_of_acidic_residues.append(
                    acid_count_this_helix)

    ttest = scipy.stats.kruskal(
        single_pass_list_of_number_of_acidic_residues, multi_pass_list_of_number_of_acidic_residues)
    n_singlepass = len(single_pass_list_of_number_of_acidic_residues)
    n_multipass = len(multi_pass_list_of_number_of_acidic_residues)
    average_singlepass = np.mean(single_pass_list_of_number_of_acidic_residues)
    average_multipass = np.mean(multi_pass_list_of_number_of_acidic_residues)

    print file, n_singlepass, n_multipass, average_singlepass, average_multipass, ttest[0], ttest[1]
