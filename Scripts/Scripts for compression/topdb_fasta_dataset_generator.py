from __future__ import division
import numpy
import sys
import re

# The input file
input_file = "top_all.txt"

# Parameters
# Due to the islice module, a try statement is used, and this means that
# only one flank length at a time can be used.
flank_lengths = [5, 10, 20]
minimum_tmd_length = 16
maximum_tmd_length = 38
flank_clash_amendment = [True, False]


for flank_length in flank_lengths:
    for flank_clash_amendment_status in flank_clash_amendment:

        output_filename = "TopDB_%s_flanklength_flankclash%s.csv" % (
            flank_length, str(flank_clash_amendment_status))
        half_flanks_output_filename = input_file.replace(
            ".txt", "_%s_flanklength_flankclash%s_only_half_flanks.csv" % (flank_length, str(flank_clash_amendment_status)))
        full_flanks_output_filename = input_file.replace(
            ".txt", "_%s_flanklength_flankclash%s_only_full_flanks.csv" % (flank_length, str(flank_clash_amendment_status)))

        number_of_records = 0
        number_of_records_correct_length = 0

        number_of_records_single = 0
        number_of_records_correct_length_single = 0
        number_of_records_multi = 0
        number_of_records_correct_length_multi = 0

        # The sequence has previously been clustered with CD-HIT to remove close
        # homologues. A list of the resulting IDs is included.
        representative_id_codes = open("ExpAll90%ID_list.txt").readlines()
        representative_id_codes = ([s.replace('\n', '')
                                    for s in representative_id_codes])

        with open(output_filename, 'w') as my_file:
            my_file.write(",TOPDB ID, N terminal inside/outside, tmh start location, tmh end location, full protein sequence, tmh sequence, N flank sequence, C flank sequence, transmembrane helix sequential number, number of transmembrane helices in protein,\n")
        my_file.closed

        short_tmds = 0
        long_tmds = 0
        tmd_counter = 0

        # The input file is parsed 3 lines at a time. The first line contains the
        # header information, the second is the fasta, and the third is the
        # topological annotation.
        records = []
        with open(input_file, "r") as database_file:
            record = []
            for line in database_file:
                record.append(line)
                if len(record) == 3:
                    records.append(record)
                    record = []
                else:
                    pass

        for this_entry in records:
            # This opens a loop that will count higher, and islice will keep
            # reading the next 3 lines until the counter exceeds the number of
            # lines available. This error is handled later and the log is outputted
            # and the script is then safely closed.

            n_terminal_start = 0

            id = this_entry[0]
            id = id.replace("\n", "")
            id = id.replace(">", "")
            sequence = this_entry[1]
            sequence = sequence.replace("\n", "")
            topology = this_entry[2]

            # This checks the ID is a representative ID.
            if id in representative_id_codes:
                id_of_record = id
                tmd_count = 0
                total_tmd_count = 0
                start_locations = []
                end_locations = []

                # This scans the topology line for the occurance when a
                # non-membrane stretch becomes a membrane (M) stretch.
                for i, c in enumerate(str(topology)):

                    if i == 0:
                        this_residue = c
                    if i > 0:
                        previous_residue = this_residue
                    this_residue = c
                    if i > 0:
                        if this_residue != previous_residue:
                            if c == "M":
                                total_tmd_count = total_tmd_count + 1
                                start_locations.append(i)
                            elif c != "M" and c != "X" and previous_residue == "M":
                                # The end locations now contains the C
                                # terminal location of each TMH.
                                end_locations.append(i)

                # Now we go through the topology and look for any flanking
                # overlaps. If a "clash" is found, the flank length is
                # adjusted accordingly.
                for i, c in enumerate(str(topology)):

                    if i == 0:
                        this_residue = c
                    if i > 0:
                        previous_residue = this_residue
                    this_residue = c

                    if i > 0:
                        if this_residue != previous_residue:
                            # The begining of TMH is recorded
                            if c == "M":

                                # This counts the number of TMDs in the
                                # sequence by looking for occurances where
                                # a non membrane annotated residue is
                                # followed by a membrane annotated residue.
                                tmd_count = tmd_count + 1

                                # Due to the process calling for a previous
                                # residue, we should start at position 1,
                                # not 0.
                                tmh_start = i + 1

                                flank1_length = flank_length
                                clash = False

                                # This counts if any of the end locations are with
                                # a flank length
                                for index in end_locations:
                                    if index < i and (index + flank_length > i - flank_length) and flank_clash_amendment_status == True:
                                        if clash == False:
                                            max_flank_size = abs(
                                                i - (index))
                                            flank1_length = max_flank_size / 2
                                            clash = True

                                # This discovers the topology of the
                                # previous residue. This allows us to
                                # assertain if the N terminal of the the
                                # TMH is inside or outside.
                                if previous_residue == "O":
                                    n_terminal_start = "Outside"
                                elif previous_residue == "I":
                                    n_terminal_start = "Inside"
                                else:
                                    pass

                            # This checks for any clashes in the C terminal flank
                            # of the tmh, or if the residue is the final residue.
                            # This also triggers the tmh to be recorded in the csv
                            # file.
                            elif (c != "M" and previous_residue == "M") or (c == "M" and i == len(str(topology))):
                                flank2_length = flank_length
                                clash = False
                                for index in reversed(start_locations):
                                    if index > i and (index - flank_length < i + flank_length) and flank_clash_amendment_status == True:
                                        if clash == False:
                                            max_flank_size = abs(
                                                i - (index))
                                            flank2_length = max_flank_size / 2
                                            clash = True

                                tmh_stop = i + 1

                                # We now have the accurate flank lengths for each
                                # tmh and the co-ordinates of the start and stop.
                                # These next statements make the strings containing
                                # the flank sequences. #-1 is to account for
                                # 0-based indexing.

                                if tmh_start - 1 - int(flank1_length) >= 0:
                                    n_terminal_flank = (
                                        sequence[tmh_start - 1 - int(flank1_length):tmh_start - 1])

                                elif tmh_start - 1 - int(flank1_length) <= 0:
                                    n_terminal_flank = (
                                        sequence[0:tmh_start - 1])

                                if tmh_stop - 1 + int(flank2_length) > len(sequence):
                                    c_terminal_flank = (
                                        sequence[tmh_stop - 1:int(len(sequence))])
                                else:
                                    c_terminal_flank = (
                                        sequence[tmh_stop - 1:(tmh_stop - 1 + int(flank2_length))])

                                tmh_sequence = sequence[
                                    tmh_start - 1:tmh_stop - 1]
                                full_sequence = sequence
                                tmd_counter = tmd_counter + 1

                                # The tmh_stop is treated with -1 to reflect
                                # how the database describes the final
                                # position. Previously in this script, it is
                                # used according to the way python slices in
                                # order to aquire the correct residue.
                                tmh_record = [id_of_record, n_terminal_start, tmh_start, tmh_stop - 1,
                                              full_sequence, tmh_sequence, n_terminal_flank, c_terminal_flank, tmd_count, total_tmd_count]

                                number_of_records = number_of_records + 1

                                if total_tmd_count == 1:
                                    number_of_records_single = number_of_records_single + 1
                                if total_tmd_count > 1:
                                    number_of_records_multi = number_of_records_multi + 1

                                if len(tmh_sequence) >= minimum_tmd_length and len(tmh_sequence) <= maximum_tmd_length:

                                    with open(output_filename, 'a') as my_file:
                                        for i in tmh_record:
                                            my_file.write(str(i))
                                            my_file.write(",")
                                        my_file.write("\n")
                                    number_of_records_correct_length = number_of_records_correct_length + 1

                                    if len(c_terminal_flank) == flank_length and len(n_terminal_flank) == flank_length:
                                        with open(full_flanks_output_filename, 'a') as my_file:
                                            my_file.write(",")
                                            for i in tmh_record:
                                                my_file.write(str(i))
                                                my_file.write(",")
                                            my_file.write("\n")

                                    if len(c_terminal_flank) >= (flank_length / 2) and len(n_terminal_flank) >= (flank_length / 2):
                                        with open(half_flanks_output_filename, 'a') as my_file:
                                            for i in tmh_record:
                                                my_file.write(str(i))
                                                my_file.write(",")
                                            my_file.write("\n")

                                    if total_tmd_count == 1:
                                        number_of_records_correct_length_single = number_of_records_correct_length_single + 1
                                    if total_tmd_count > 1:
                                        number_of_records_correct_length_multi = number_of_records_correct_length_multi + 1

        # Once there are no more lines in the input file, the following logs are
        # printed.
        this_entry = 'null'
        print "Records"
        print "Total records", number_of_records
        print "Total records after length filter", number_of_records_correct_length
        print "Single-pass"
        print "Total:", number_of_records_single
        print "...after length filter:", number_of_records_correct_length_single
        print "Multi-pass"
        print "Total:", number_of_records_multi
        print "...after length filter:", number_of_records_correct_length_multi
