from __future__ import division


list_of_files = [
    "TopDB_20_flanklength_flankclashTrue.csv",
    "UniHuman_20_flanklength_flankclashTrue.csv",
    ]


for file in list_of_files:
    charge_multipass = []
    charge_singlepass = []

    results = []

    with open(file) as inputfile:
        for line in inputfile:
            results.append(line.strip().split(','))

    maximum_tmd_length = 0
    for entry in results:
        if entry == results[0]:
            pass
        else:
            tmh_sequence = entry[6]
            if len(tmh_sequence) > maximum_tmd_length:
                maximum_tmd_length = len(tmh_sequence)

    max_sequence_length = maximum_tmd_length + 40

    for _ in range(max_sequence_length + 1):
        charge_multipass.append(0)
        charge_singlepass.append(0)

    list_of_tmh_segments = []
    list_of_singlepass_tmh_segments = []
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
                tmh_unaltered_sequence = str(
                    n_flank_sequence) + str(tmh_sequence) + str(c_flank_sequence)
                tmh_reversed_sequence = tmh_unaltered_sequence[::-1]
                correction_number = (
                    (maximum_tmd_length - len(tmh_sequence)) / 2 + (20 - len(c_flank_sequence)))
                tmh_segment = "J" * int(correction_number)
                tmh_segment = tmh_segment + tmh_reversed_sequence
            if "Inside" in str(n_terminal_start):
                tmh_unaltered_sequence = str(
                    n_flank_sequence) + str(tmh_sequence) + str(c_flank_sequence)
                correction_number = (
                    (maximum_tmd_length - len(tmh_sequence)) / 2 + (20 - len(n_flank_sequence)))
                tmh_segment = "J" * int(correction_number)
                tmh_segment = tmh_segment + tmh_unaltered_sequence
            list_of_tmh_segments.append(tmh_segment)

            if int(total_tmd_count) == 1:
                list_of_singlepass_tmh_segments.append(tmh_segment)
            elif int(total_tmd_count) > 1:
                list_of_multipass_tmh_segments.append(tmh_segment)

    for i in list_of_multipass_tmh_segments:
        for position, residue in enumerate(str(i)):
            if residue == "R" or residue == "K":
                charge_multipass[position] = charge_multipass[position] + 1
            if residue == "E" or residue == "D":
                charge_multipass[position] = charge_multipass[position] - 1

    for i in list_of_singlepass_tmh_segments:
        for position, residue in enumerate(str(i)):
            if residue == "R" or residue == "K":
                charge_singlepass[position] = charge_singlepass[position] + 1
            if residue == "E" or residue == "D":
                charge_singlepass[position] = charge_singlepass[position] - 1

    for position, charge in enumerate(charge_singlepass):
        relative_charge = (charge / len(list_of_singlepass_tmh_segments))
        charge_singlepass[position] = relative_charge

    for position, charge in enumerate(charge_multipass):
        relative_charge = (charge / len(list_of_multipass_tmh_segments))
        charge_multipass[position] = relative_charge

    print file
    '''print "Single-pass:", len(list_of_singlepass_tmh_segments)
    print "Multi-pass:", len(list_of_multipass_tmh_segments)
    print "Combined:", len(list_of_tmh_segments)'''
    print "Single-pass charge profile:", charge_singlepass
    print "Multi-pass charge profile:", charge_multipass
