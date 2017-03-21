from __future__ import division

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


print "file, single_inside, single_outside, multi_inside, multi_outside"

for file in list_of_files:
    single_inside = 0
    single_outside = 0
    multi_inside = 0
    multi_outside = 0

    results = []
    with open(file) as inputfile:
        for line in inputfile:
            results.append(line.strip().split(','))
    previous_id="NULL"
    for entry in results:
        if entry == results[0]:
            pass
        else:
            name = entry[0]
            id = entry[1]
            n_terminal_start = str(entry[2])
            tmh_start_location = entry[3]
            tmh_end_location = entry[4]
            sequence = entry[5]
            tmh_sequence = entry[6]
            N_flank_sequence = entry[7]
            C_flank_sequence = entry[8]




            tmh_number = int(entry[9])
            total_tmd_count = int(entry[10])

            if tmh_number == 1:
                if total_tmd_count == 1:
                    if n_terminal_start == "Inside":
                        single_inside = single_inside + 1
                    elif n_terminal_start == "Outside":
                        single_outside = single_outside + 1
                    else:
                        print "Error. No IO found."
                elif total_tmd_count > 1:
                    if n_terminal_start == "Inside":
                        multi_inside = multi_inside + 1
                    elif n_terminal_start == "Outside":
                        multi_outside = multi_outside + 1
                    else:
                        print "Error. No IO found."
                else:
                    print "Error. Incompatible tmd count."
            elif previous_id != id:
                #This means that the "first" id recorded in a record is no the first TMH.
                #print id
                pass
            else:
                pass
            previous_id = id
    print file,",", single_inside,",",(single_inside/(single_inside+single_outside))*100,",", single_outside, ",",(single_outside/(single_inside+single_outside))*100,",", multi_inside, ",",(multi_inside/(multi_inside+multi_outside))*100,",", multi_outside, ",",(multi_outside/(multi_inside+multi_outside))*100
