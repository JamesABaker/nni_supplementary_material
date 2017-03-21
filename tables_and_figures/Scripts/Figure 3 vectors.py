from __future__ import division
import numpy as np
import pylab
import itertools


# Add files to convert from csv
list_of_files = ["UniHuman_20_flanklength.csv",
"TopDB_20_flanklength.csv"]




residues = ["I", "V", "L", "F", "C", "M", "A", "G", "T", "S", "W", "Y", "P", "H", "E", "Q", "D", "N", "K", "R"]

for file in list_of_files:
    list_to_print_singlepass = []
    list_to_print_multipass = []
    # results variable contains a list of lists of each line in the csv as a record.
    results=[]
    with open(file) as inputfile:
        for line in inputfile:
            results.append(line.strip().split(','))


    #The total list of residues is held in a list. This will be turned into a string for counting later.
    list_of_slices=[]

    # Values for longest TMD
    maximum_tmd_length = 0
    maximum_segment_length =0
    maxCflank_length=0
    maxNflank_length=0
    maxflank_length=0
    that_length = 0

    #Values for longest inside half for alignment to point 0
    maximum_tmd_length_for_alignment = 0
    maximum_inner_flank_length = 0

    for entry in results:
        # Don't read header line
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

            #Finding where to align 0 to.
            #Use the inner flank length and the inner tmh leaflet (1/2 TMD) and find the longest of that combo.
            if "Outside" in str(n_terminal_start):
                inner_flank_length = len(C_flank_sequence)
            if "Inside" in str(n_terminal_start):
                inner_flank_length = len(N_flank_sequence)
            tmd_length=len(tmh_sequence)

            if tmd_length > maximum_tmd_length_for_alignment:
                maximum_tmd_length_for_alignment = tmd_length
            if inner_flank_length > maximum_inner_flank_length:
                maximum_inner_flank_length = inner_flank_length

            #Finding the longest sequence:
            this_tmh_and_flank_length = len(tmh_sequence)+len(N_flank_sequence)+len(C_flank_sequence)
            if this_tmh_and_flank_length > maximum_segment_length:
                maximum_length = len(tmh_sequence)+len(N_flank_sequence)+len(C_flank_sequence)
                maximum_tmd_length = len(tmh_sequence)

                if "Outside" in str(n_terminal_start):
                    tmh_unaltered_sequence = str(N_flank_sequence) + str(tmh_sequence) + str(C_flank_sequence)
                    tmh_reversed_sequence = tmh_unaltered_sequence[::-1]
                    maxflank_length = len(C_flank_sequence)

                if "Inside" in str(n_terminal_start):
                    tmh_unaltered_sequence = str(N_flank_sequence) + str(tmh_sequence) + str(C_flank_sequence)
                    maxflank_length = len(N_flank_sequence)

                maximum_segment_length = this_tmh_and_flank_length

    max_sequence_length = maximum_segment_length

    list_of_tmh_segments = []
    list_of_singlepass_tmh_segments = []
    list_of_multipass_tmh_segments = []

    for entry in results:
        # This is the column title line
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

            #Aligning the residues according to the reference TMH (longest inner flank and inner TMH). Effictively, they are all aligned to 0.
            if "Outside" in str(n_terminal_start):
                tmh_unaltered_sequence = str(N_flank_sequence) + str(tmh_sequence) + str(C_flank_sequence)
                tmh_reversed_sequence = tmh_unaltered_sequence[::-1]
                correction_number = ((maximum_tmd_length_for_alignment/2) + maximum_inner_flank_length) - ((len(tmh_sequence)/ 2) + len(C_flank_sequence))
                tmh_segment = "J" * (int(correction_number) + 1) # +1 used since python is counting from 0, however I'm not sure entirely that it's needed here.
                tmh_segment = tmh_segment + tmh_reversed_sequence
            if "Inside" in str(n_terminal_start):
                tmh_unaltered_sequence = str(N_flank_sequence) + str(tmh_sequence) + str(C_flank_sequence)
                correction_number = ((maximum_tmd_length_for_alignment/2) + maximum_inner_flank_length) - ((len(tmh_sequence)/ 2) + len(N_flank_sequence))
                tmh_segment = "J" * (int(correction_number) + 1)
                tmh_segment = tmh_segment + tmh_unaltered_sequence
            list_of_tmh_segments.append(tmh_segment)

            # Builds up list of sequences in order with dummy "Js" as place
            # holders (these are ignored in counting, however act as
            # placeholders so the sequences can be counted)
            if int(total_tmd_count) == 1:
                list_of_singlepass_tmh_segments.append(tmh_segment)
            if int(total_tmd_count) > 1:
                list_of_multipass_tmh_segments.append(tmh_segment)

    for residue in residues:

        residue_singlepass = [0]
        residue_multipass = [0]

        for _ in range(max_sequence_length):
            residue_singlepass.append(0)
            residue_multipass.append(0)

        # Counting residues at each position for singlepass.
        for i in list_of_singlepass_tmh_segments:
            for position, each_residue in enumerate(str(i)):
                if each_residue == residue:
                    try:
                        residue_singlepass[position] = residue_singlepass[position] + 1
                    except IndexError:
                        pass

        # Counting residues at each position for multipass.
        for i in list_of_multipass_tmh_segments:
            for position, each_residue in enumerate(str(i)):
                if each_residue == residue:
                    try:
                        residue_multipass[position] = residue_multipass[position] + 1
                    except IndexError:
                        pass


        #Finding the normalised values
        residue_singlepass_normalised = []
        for i in residue_singlepass:
            try:
                normalised_value = (100 * i)/ sum(residue_singlepass)
                residue_singlepass_normalised.append(normalised_value)
            except ZeroDivisionError:
                residue_singlepass_normalised.append(0)
        tmh_vectors_to_print = [list(residue), list(residue_singlepass_normalised)]
        list_to_print_singlepass.append(tmh_vectors_to_print)

        #Finding the normalised values
        residue_multipass_normalised = []
        for i in residue_multipass:
            try:
                normalised_value = (100 * i)/ sum(residue_multipass)
                residue_multipass_normalised.append(normalised_value)
            except ZeroDivisionError:
                residue_multipass_normalised.append(0)
        tmh_vectors_to_print = [list(residue), list(residue_multipass_normalised)]
        list_to_print_multipass.append(tmh_vectors_to_print)


    #sequence position for the x axis.
    sequence_position=[]
    for position_number, item in enumerate(residue_multipass):
        sequence_number=position_number -(max_sequence_length/2)
        sequence_position.append(sequence_number)

    print "\n", file
    print "\nSequence position,", sequence_position
    #output

    print "\nSinglepass"
    for vector in list_to_print_singlepass:
        print vector

    print "\nMultipass"
    for vector in list_to_print_multipass:
        print vector
