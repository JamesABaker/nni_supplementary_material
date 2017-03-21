from __future__ import division
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pylab
import itertools
from time import gmtime, strftime
from random import random, seed
from matplotlib import cm
from matplotlib import rcParams



list_of_files = [#"TopDB_10_flanklength_flankclashFalse.csv",
                 "TopDB_10_flanklength_flankclashTrue.csv",
                 #"UniHuman_10_flanklength_flankclashFalse.csv",
                 "UniHuman_10_flanklength_flankclashTrue.csv",
                 ]


with open('List_of_Complex_ids.txt') as f:
    complex_id_list = f.read().splitlines()

with open('List_of_Simple_ids.txt') as f:
    simple_id_list = f.read().splitlines()

list_of_residues = [["C"], ["R"], ["K"], ["R","K"], ["D"], ["E"], ["D","E"], ["Y"], ["W"], ["L"]]

for file in list_of_files:
    print file
    for residue_type in list_of_residues:

        results=[]
        with open(file) as inputfile:
            for line in inputfile:
                results.append(line.strip().split(','))

        list_of_slices=[]

        singlepass_helices = 0

        residue_simple = [0]
        residue_complex = [0]
        residue_multipass =[0]

        maximum_tmd_length = 0
        maximum_segment_length =0
        maxflank_length=0
        that_length = 0

        maximum_tmd_length_for_alignment = 0
        maximum_inner_flank_length = 0

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

                if "Outside" in str(n_terminal_start):
                    inner_flank_length = len(c_flank_sequence)
                if "Inside" in str(n_terminal_start):
                    inner_flank_length = len(n_flank_sequence)
                tmd_length=len(tmh_sequence)

                if tmd_length > maximum_tmd_length_for_alignment:
                    maximum_tmd_length_for_alignment = tmd_length
                if inner_flank_length > maximum_inner_flank_length:
                    maximum_inner_flank_length = inner_flank_length

                this_tmh_and_flank_length = len(tmh_sequence)+len(n_flank_sequence)+len(c_flank_sequence)
                if this_tmh_and_flank_length > maximum_segment_length:
                    maximum_length = len(tmh_sequence)+len(n_flank_sequence)+len(c_flank_sequence)
                    maximum_tmd_length = len(tmh_sequence)

                    if "Outside" in str(n_terminal_start):
                        tmh_unaltered_sequence = str(n_flank_sequence) + str(tmh_sequence) + str(c_flank_sequence)
                        tmh_reversed_sequence = tmh_unaltered_sequence[::-1]
                        maxflank_length = len(c_flank_sequence)

                    if "Inside" in str(n_terminal_start):
                        tmh_unaltered_sequence = str(n_flank_sequence) + str(tmh_sequence) + str(c_flank_sequence)
                        maxflank_length = len(n_flank_sequence)

                    maximum_segment_length = this_tmh_and_flank_length

        max_sequence_length = maximum_segment_length

        for _ in range(max_sequence_length):
            residue_simple.append(0)
            residue_complex.append(0)
            residue_multipass.append(0)

        list_of_tmh_segments = []
        list_of_simple_singlepass_tmh_segments = []
        list_of_complex_singlepass_tmh_segments = []
        list_of_multipass_tmh_segments = []

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
                    tmh_unaltered_sequence = str(n_flank_sequence) + str(tmh_sequence) + str(c_flank_sequence)
                    tmh_reversed_sequence = tmh_unaltered_sequence[::-1]
                    correction_number = ((maximum_tmd_length_for_alignment/2) + maximum_inner_flank_length) - ((len(tmh_sequence)/ 2) + len(c_flank_sequence))
                    tmh_segment = "J" * (int(correction_number) + 1)
                    tmh_segment = tmh_segment + tmh_reversed_sequence

                if "Inside" in str(n_terminal_start):
                    tmh_unaltered_sequence = str(n_flank_sequence) + str(tmh_sequence) + str(c_flank_sequence)
                    correction_number = ((maximum_tmd_length_for_alignment/2) + maximum_inner_flank_length) - ((len(tmh_sequence)/ 2) + len(n_flank_sequence))
                    tmh_segment = "J" * (int(correction_number) + 1)
                    tmh_segment = tmh_segment + tmh_unaltered_sequence
                list_of_tmh_segments.append(tmh_segment)

                if int(total_tmd_count) == 1 and str(id) in simple_id_list:
                    list_of_simple_singlepass_tmh_segments.append(tmh_segment)

                if int(total_tmd_count) == 1 and str(id) in complex_id_list:
                    list_of_complex_singlepass_tmh_segments.append(tmh_segment)

                if int(total_tmd_count) > 1:
                    list_of_multipass_tmh_segments.append(tmh_segment)


        print "Simple single-pass,", len(list_of_simple_singlepass_tmh_segments)
        print "Complex single-pass,", len(list_of_complex_singlepass_tmh_segments)

        for i in list_of_simple_singlepass_tmh_segments:
            for position, residue in enumerate(str(i)):
                if residue in residue_type:
                    try:
                        residue_simple[position] = residue_simple[position] + 1
                    except IndexError:
                        pass
                if residue == "J":
                    pass

        for i in list_of_complex_singlepass_tmh_segments:
            for position, residue in enumerate(str(i)):
                if residue in residue_type:
                    try:
                        residue_complex[position] = residue_complex[position] + 1
                    except IndexError:
                        pass
                if residue == "J":
                    pass

        for i in list_of_multipass_tmh_segments:
            for position, residue in enumerate(str(i)):
                if residue in residue_type:
                    try:
                        residue_multipass[position] = residue_multipass[position] + 1
                    except IndexError:
                        pass
                if residue == "J":
                    pass

        list(itertools.chain.from_iterable(list_of_simple_singlepass_tmh_segments))

        residue_simple_normalised = []
        for i in residue_simple:
            try:
                normalised_value = (100 * i)/ sum(residue_simple)
                residue_simple_normalised.append(normalised_value)
            except ZeroDivisionError:
                residue_simple_normalised.append(0)

        residue_complex_normalised = []
        for i in residue_complex:
            try:
                normalised_value = (100 * i)/ sum(residue_complex)
                residue_complex_normalised.append(normalised_value)
            except ZeroDivisionError:
                residue_complex_normalised.append(0)

        residue_multipass_normalised = []
        for i in residue_multipass:
            try:
                normalised_value = (100 * i)/ sum(residue_multipass)
                residue_multipass_normalised.append(normalised_value)
            except ZeroDivisionError:
                residue_multipass_normalised.append(0)



        sequence_position=[]
        for position_number, item in enumerate(residue_multipass):
            sequence_number=position_number -(max_sequence_length/2)
            sequence_position.append(sequence_number)

        print file
        print residue_type
        print sequence_position

        print "Simple,", residue_simple
        print "Complex,",residue_complex
        print "Multi-pass,", residue_multipass
        print "\n"

        ax = plt.subplot(111)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

        ax.add_patch(patches.Rectangle(
                (-10, 0),   # (x,y)
                20,          # width
                100,          # height
                alpha=0.01
            )
        )
        ax.add_patch(patches.Rectangle(
                (-30, 0),   # (x,y)
                5,          # width
                100,          # height
                alpha=0.1,
                facecolor="gray"
            )
        )
        ax.add_patch(patches.Rectangle(
                (25, 0),   # (x,y)
                5,          # width
                100,          # height
                alpha=0.1,
                facecolor="gray"
            )
        )

        ax.add_patch(patches.Rectangle(
                (-20, 0),   # (x,y)
                10,          # width
                100,          # height
                alpha=0.1,
                facecolor="gray"
            )
        )
        ax.add_patch(patches.Rectangle(
                (10, 0),   # (x,y)
                10,          # width
                100,          # height
                alpha=0.1,
                facecolor="gray"
            )
        )
        plt.plot((0, 0), (0, 100), linestyle='-', linewidth=1,color='black')
        plt.plot(sequence_position, residue_simple_normalised, linestyle='-', marker='',  linewidth=2,  color='lightskyblue')
        plt.plot(sequence_position, residue_complex_normalised, linestyle='-', marker='',  linewidth=2,  color='lightcoral')
        plt.plot(sequence_position, residue_multipass_normalised, linestyle='-', marker='',  linewidth=2,  color='black')


        font = {'size':24}

        if residue_type == ["D","E"]:
            residue_title = "D and E"
        elif residue_type == ["R","K"]:
            residue_title = "R and K"
        else:
            residue_title = residue_type[0]

        plt.xlabel('Sequence Position', **font)
        plt.ylabel('%s Relative Percentage' % residue_title, **font)
        plt.tick_params(labelsize=24)
        pylab.xlim([-20,20])
        pylab.ylim([0, 10])
        timestamp = str(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        filename = timestamp + file.replace('.csv', '') + residue_title + ".pdf"
        plt.gcf().subplots_adjust(bottom=0.2)
        plt.savefig(filename)
        #plt.show()

        plt.clf()
        plt.cla()
