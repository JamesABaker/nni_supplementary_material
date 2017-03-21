from __future__ import division
import subprocess
import re
import scipy
from scipy import stats
import numpy as np
import pylab as P
import random
import ast


header = ">Token_header"



list_of_files = ["UniER_20_flanklength_flankclashTrue.csv"]


list_of_scales = ["KyteDoolittle.pl"]

for scale in list_of_scales:
    for file in list_of_files:
        pvaluelist = []

        results = []

        with open(file) as inputfile:
            for line in inputfile:
                results.append(line.strip().split(','))

        maximum_tmd_length = 0
        print "Identifying longest sequence"

        number_of_entries = 0
        for entry in results:
            number_of_entries = number_of_entries + 1

            if entry == results[0]:
                pass
            else:
                tmh_sequence = entry[6]
                if len(tmh_sequence) > maximum_tmd_length:
                    maximum_tmd_length = len(tmh_sequence)

        max_sequence_length = maximum_tmd_length + 40

        single = []
        multi = []
        entries_calculated_so_far = 0

        print "Parsing csv file and calculating windowed hydrophobicity..."
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
                    tmh_unaltered_sequence = str(
                        n_flank_sequence) + str(tmh_sequence) + str(c_flank_sequence)
                    tmh_reversed_sequence = tmh_unaltered_sequence[::-1]
                    correction_number = (
                        (maximum_tmd_length - len(tmh_sequence)) / 2 + (20 - len(c_flank_sequence)))
                    tmh_segment = tmh_reversed_sequence

                if "Inside" in str(n_terminal_start):
                    tmh_unaltered_sequence = str(
                        n_flank_sequence) + str(tmh_sequence) + str(c_flank_sequence)
                    correction_number = (
                        (maximum_tmd_length - len(tmh_sequence)) / 2 + (20 - len(n_flank_sequence)))
                    tmh_segment = tmh_unaltered_sequence

                sequence = tmh_segment

                with open('KD_calc_in.txt', 'wb') as temp_fasta:
                    temp_fasta.write(header)
                    temp_fasta.write("\n")
                    temp_fasta.write(sequence)

                var = "/"
                pipe = subprocess.Popen(["perl", scale, var])
                pipe.wait()

                with open('KDcalc_out.txt', 'rb') as totalKD:
                    lines = totalKD.readlines()
                    KD_line = lines[4]
                    totalKD.close()

                    KD = str(KD_line)

                    KD = re.sub("[KD=]", '', KD)
                    result = KD
                    result = re.sub("  ", ', ', result)
                    result = re.sub(", , , ,", '', result)
                    result = re.sub(" , ", '', result)
                    result = re.sub("\n", '', result)

                result = result.replace(" -", ", -")
                result = result.replace(",, -", ", -")
                hydrophobicity = "," * int(correction_number)
                hydrophobicity = hydrophobicity + result
                output_line = [id, total_tmd_count,
                               correction_number, hydrophobicity]

                x = "[" + str(output_line[3]) + "]"
                x = x.replace(",,", ",'',")
                x = x.replace(",,", ",'',")
                x = x.replace("[,", "['',")

                values = ast.literal_eval(x)

                while len(values) < max_sequence_length:
                    values.append("")
                if int(total_tmd_count) == 1:
                    print "single appended", total_tmd_count
                    single.append(values)
                if int(total_tmd_count) > 1:
                    print "multi appended", total_tmd_count
                    multi.append(values)

                entries_calculated_so_far = entries_calculated_so_far + 1
                print entries_calculated_so_far, "/", number_of_entries

        print "Calculating pvalues."


        for i in range(int(max_sequence_length)-1):
            values_for_single = []
            values_for_multi = []

            for n in single:
                if n[i] == "":
                    pass
                else:
                    values_for_single.append(float(n[i]))

            #print values_for_single
            for y in multi:
                #integer_check = type( y[i] ) == int
                if y[i] == "":
                    pass
                else:
                    values_for_multi.append(float(y[i]))

            #print values_for_multi
            if len(values_for_multi) > 0 and len(values_for_single) > 0:
                p_value_at_this_position = scipy.stats.kruskal(values_for_single, values_for_multi)
                pvaluelist.append(p_value_at_this_position[1])
            else:
                print "No values at position. Using '1' as a placeholder."
                pvaluelist.append(1)

        print scale, ",", pvaluelist
