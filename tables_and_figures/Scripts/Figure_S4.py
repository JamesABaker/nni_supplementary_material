from __future__ import division
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pylab
import itertools
from scipy import stats
from time import gmtime, strftime
from random import random, seed
from matplotlib import cm
from matplotlib import rcParams
from matplotlib.ticker import FuncFormatter

matplotlib.rcParams.update({'font.size': 22})


list_of_files = [
    "TopDB_10_flanklength_flankclashFalse.csv",
    "TopDB_10_flanklength_flankclashTrue.csv",
    "TopDB_20_flanklength_flankclashFalse.csv",
    "TopDB_20_flanklength_flankclashTrue.csv",
    "TopDB_5_flanklength_flankclashFalse.csv",
    "TopDB_5_flanklength_flankclashTrue.csv",
    "UniArch_10_flanklength_flankclashFalse.csv",
    "UniArch_10_flanklength_flankclashTrue.csv",
    "UniArch_20_flanklength_flankclashFalse.csv",
    "UniArch_20_flanklength_flankclashTrue.csv",
    "UniArch_5_flanklength_flankclashFalse.csv",
    "UniArch_5_flanklength_flankclashTrue.csv",
    "UniBacilli_10_flanklength_flankclashFalse.csv",
    "UniBacilli_10_flanklength_flankclashTrue.csv",
    "UniBacilli_20_flanklength_flankclashFalse.csv",
    "UniBacilli_20_flanklength_flankclashTrue.csv",
    "UniBacilli_5_flanklength_flankclashFalse.csv",
    "UniBacilli_5_flanklength_flankclashTrue.csv",
    "UniCress_10_flanklength_flankclashFalse.csv",
    "UniCress_10_flanklength_flankclashTrue.csv",
    "UniCress_20_flanklength_flankclashFalse.csv",
    "UniCress_20_flanklength_flankclashTrue.csv",
    "UniCress_5_flanklength_flankclashFalse.csv",
    "UniCress_5_flanklength_flankclashTrue.csv",
    "UniEcoli_10_flanklength_flankclashFalse.csv",
    "UniEcoli_10_flanklength_flankclashTrue.csv",
    "UniEcoli_20_flanklength_flankclashFalse.csv",
    "UniEcoli_20_flanklength_flankclashTrue.csv",
    "UniEcoli_5_flanklength_flankclashFalse.csv",
    "UniEcoli_5_flanklength_flankclashTrue.csv",
    "UniER_10_flanklength_flankclashFalse.csv",
    "UniER_10_flanklength_flankclashTrue.csv",
    "UniER_20_flanklength_flankclashFalse.csv",
    "UniER_20_flanklength_flankclashTrue.csv",
    "UniER_5_flanklength_flankclashFalse.csv",
    "UniER_5_flanklength_flankclashTrue.csv",
    "UniFungi_10_flanklength_flankclashFalse.csv",
    "UniFungi_10_flanklength_flankclashTrue.csv",
    "UniFungi_20_flanklength_flankclashFalse.csv",
    "UniFungi_20_flanklength_flankclashTrue.csv",
    "UniFungi_5_flanklength_flankclashFalse.csv",
    "UniFungi_5_flanklength_flankclashTrue.csv",
    "UniGolgi_10_flanklength_flankclashFalse.csv",
    "UniGolgi_10_flanklength_flankclashTrue.csv",
    "UniGolgi_20_flanklength_flankclashFalse.csv",
    "UniGolgi_20_flanklength_flankclashTrue.csv",
    "UniGolgi_5_flanklength_flankclashFalse.csv",
    "UniGolgi_5_flanklength_flankclashTrue.csv",
    "UniHuman_10_flanklength_flankclashFalse.csv",
    "UniHuman_10_flanklength_flankclashTrue.csv",
    "UniHuman_20_flanklength_flankclashFalse.csv",
    "UniHuman_20_flanklength_flankclashTrue.csv",
    "UniHuman_5_flanklength_flankclashFalse.csv",
    "UniHuman_5_flanklength_flankclashTrue.csv",
    "UniPM_10_flanklength_flankclashFalse.csv",
    "UniPM_10_flanklength_flankclashTrue.csv",
    "UniPM_20_flanklength_flankclashFalse.csv",
    "UniPM_20_flanklength_flankclashTrue.csv",
    "UniPM_5_flanklength_flankclashFalse.csv",
    "UniPM_5_flanklength_flankclashTrue.csv",
]
list_of_files = [
    "TopDB_20_flanklength_flankclashFalse.csv",
    "TopDB_20_flanklength_flankclashTrue.csv",
    "UniHuman_20_flanklength_flankclashFalse.csv",
    "UniHuman_20_flanklength_flankclashTrue.csv",
]


type_of_proteins = [
    "singlepass",
    "multipass",
]


def single_or_multi(type_of_protein, total_tmd_count, record_id):
    if type_of_protein == "singlepass" and total_tmd_count == 1:
        return True
    elif type_of_protein == "multipass" and total_tmd_count > 1:
        return True
    else:
        return False


for file in list_of_files:
    for type_of_protein in type_of_proteins:
        tmh_count = []
        inside_flank_count = []
        outside_flank_count = []

        results = []
        with open(file) as inputfile:
            for line in inputfile:
                results.append(line.strip().split(','))

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

                if single_or_multi(type_of_protein, int(total_tmd_count), id) == True:
                    if n_terminal_start == "Outside":
                        inside_length = len(c_flank_sequence)
                        tmh_length = len(tmh_sequence)
                        outside_length = len(n_flank_sequence)

                    elif n_terminal_start == "Inside":
                        inside_length = len(n_flank_sequence)
                        tmh_length = len(tmh_sequence)
                        outside_length = len(c_flank_sequence)

                    else:
                        pass

                    tmh_count.append(tmh_length)
                    inside_flank_count.append(inside_length)
                    outside_flank_count.append(outside_length)

        features_to_plot = [inside_flank_count, tmh_count, outside_flank_count]

        for list_to_plot in features_to_plot:

            if list_to_plot == tmh_count:
                feature_type = "TMHs"
            elif list_to_plot == inside_flank_count:
                feature_type = "inside flanks"
            elif list_to_plot == outside_flank_count:
                feature_type = "outside flanks"
            else:
                pass

            timestamp = str(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
            filename = timestamp + \
                file.replace('.csv', '') + "_" + type_of_protein + \
                "_" + feature_type + ".pdf"

            def to_percent(y, position):
                '''
                Scales the axis to a percentage.
                '''
                s = str(int(y / len(list_to_plot) * 100))

                # The percent symbol needs escaping in latex
                if matplotlib.rcParams['text.usetex'] is True:
                    return s + r'$\%$'
                else:
                    return s + '%'
            a = np.array(list_to_plot)
            counts = np.bincount(a)
            mode = np.argmax(counts)
            counter = 0
            for i in list_to_plot:
                if mode == i:
                    counter = counter + 1

            print "Making", filename
            plt.hist(list_to_plot)
            plt.gca().set_ylim([0, len(list_to_plot)])
            plt.yticks(np.arange(0, len(list_to_plot), len(list_to_plot) / 4))
            formatter = FuncFormatter(to_percent)
            plt.gca().yaxis.set_major_formatter(formatter)
            plt.gca().set_xlim([int(min(list_to_plot)), int(max(list_to_plot))])
            plt.xlabel('Length in residues of %s' % feature_type)
            plt.ylabel('Percentage')
            plt.tight_layout()
            plt.savefig(filename)
            # plt.show()
            plt.clf()
            plt.cla()
