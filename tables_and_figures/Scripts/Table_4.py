from __future__ import division
import scipy
from scipy import stats
import numpy as np
import pylab as P
import random

print "Leucine residue statistics.\n"
print "File, SP leucine residues inside, SP leucine residues outside, SP Percentage difference, SP t-statistic, SP P-value, MP negative leucine inside, MP leucine residues outside, MP percentage difference, MP t-statstic, MP P-value"

list_of_files = ["TopDB_10_flanklength_flankclashFalse.csv",
                 "TopDB_10_flanklength_flankclashTrue.csv",
                 "UniHuman_10_flanklength_flankclashFalse.csv",
                 "UniHuman_10_flanklength_flankclashTrue.csv",
                 "UniER_10_flanklength_flankclashFalse.csv",
                 "UniER_10_flanklength_flankclashTrue.csv",
                 "UniGolgi_10_flanklength_flankclashFalse.csv",
                 "UniGolgi_10_flanklength_flankclashTrue.csv",
                 "UniPM_10_flanklength_flankclashFalse.csv",
                 "UniPM_10_flanklength_flankclashTrue.csv",
                 "UniCress_10_flanklength_flankclashFalse.csv",
                 "UniCress_10_flanklength_flankclashTrue.csv",
                 "UniFungi_10_flanklength_flankclashFalse.csv",
                 "UniFungi_10_flanklength_flankclashTrue.csv",
                 "UniBacilli_10_flanklength_flankclashFalse.csv",
                 "UniBacilli_10_flanklength_flankclashTrue.csv",
                 "UniEcoli_10_flanklength_flankclashFalse.csv",
                 "UniEcoli_10_flanklength_flankclashTrue.csv",
                 "UniArch_10_flanklength_flankclashFalse.csv",
                 "UniArch_10_flanklength_flankclashTrue.csv",
                 ]


for file in list_of_files:

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
            n_terminal_start = entry[2]
            tmh_start_location = entry[3]
            tmh_end_location = entry[4]
            sequence = entry[5]
            tmh_sequence = entry[6]
            n_flank_sequence = entry[7]
            c_flank_sequence = entry[8]
            tmh_number = entry[9]
            total_tmd_count = int(entry[10])
            correction_number = 0

            tmh_length=len(tmh_sequence)
            N_leaflet=tmh_sequence[0:int((tmh_length/2))]
            C_leaflet=tmh_sequence[-int((tmh_length/2)):]

            if "Outside" in str(n_terminal_start):
                inside_leaflet = C_leaflet[::-1]
                outside_leaflet = N_leaflet[::-1]

            if "Inside" in str(n_terminal_start):
                inside_leaflet = N_leaflet
                outside_leaflet = C_leaflet



            if int(total_tmd_count) == 1:
                inside_leucine_count = inside_leaflet.count('L')
                outside_leucine_count = outside_leaflet.count('L')

                single_inside_leaflets.append(inside_leucine_count)
                single_outside_leaflets.append(outside_leucine_count)

                single_total_leucines_leaflets = single_total_leucines_leaflets + inside_leucine_count +outside_leucine_count


            elif int(total_tmd_count) > 1:
                inside_leucine_count = inside_leaflet.count('L')
                outside_leucine_count = outside_leaflet.count('L')

                multi_inside_leaflets.append(inside_leucine_count)
                multi_outside_leaflets.append(outside_leucine_count)

                multi_total_leucines_leaflets = multi_total_leucines_leaflets + inside_leucine_count + outside_leucine_count

    ttest_singlepass= scipy.stats.kruskal(single_inside_leaflets, single_outside_leaflets)
    percentage_singlepass = np.sum(single_inside_leaflets)/np.sum(single_outside_leaflets)*100
    percentage_multipass = np.sum(multi_inside_leaflets)/np.sum(multi_outside_leaflets)*100

    ttest_multipass= scipy.stats.kruskal(multi_inside_leaflets, multi_outside_leaflets)

    print file, np.sum(single_inside_leaflets), np.sum(single_outside_leaflets),percentage_singlepass, ttest_singlepass[0], ttest_singlepass[1], np.sum(multi_inside_leaflets), percentage_multipass,  np.sum(multi_outside_leaflets), ttest_multipass[0], ttest_multipass[1]
