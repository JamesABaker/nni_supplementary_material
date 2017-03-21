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


list_of_files = ["TopDB_10_flanklength_flankclashFalse.csv",
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

for file in list_of_files:
    print "\n",file
    results=[]
    with open(file) as inputfile:
        for line in inputfile:
            results.append(line.strip().split(','))

    total_singlepass_inner_insidefraction_sequences = []
    total_singlepass_outer_insidefraction_sequences = []

    total_singlepass_inner_outsidefraction_sequences = []
    total_singlepass_outer_outsidefraction_sequences = []

    list_of_slices=[]

    singlepass_helices = 0

    A_singlepass = [0]
    R_singlepass = [0]
    N_singlepass = [0]
    D_singlepass = [0]
    B_singlepass = [0]
    C_singlepass = [0]
    E_singlepass = [0]
    Q_singlepass = [0]
    Z_singlepass = [0]
    G_singlepass = [0]
    H_singlepass = [0]
    I_singlepass = [0]
    L_singlepass = [0]
    K_singlepass = [0]
    M_singlepass = [0]
    F_singlepass = [0]
    P_singlepass = [0]
    S_singlepass = [0]
    T_singlepass = [0]
    W_singlepass = [0]
    Y_singlepass = [0]
    V_singlepass = [0]

    maximum_tmd_length = 0
    maximum_segment_length =0
    maxCflank_length=0
    maxNflank_length=0
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
            N_flank_sequence = entry[7]
            C_flank_sequence = entry[8]
            tmh_number = entry[9]
            total_tmd_count = entry[10]

            if "Outside" in str(n_terminal_start):
                inner_flank_length = len(C_flank_sequence)
            if "Inside" in str(n_terminal_start):
                inner_flank_length = len(N_flank_sequence)
            tmd_length=len(tmh_sequence)

            if tmd_length > maximum_tmd_length_for_alignment:
                maximum_tmd_length_for_alignment = tmd_length
            if inner_flank_length > maximum_inner_flank_length:
                maximum_inner_flank_length = inner_flank_length

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

    for _ in range(max_sequence_length):
        A_singlepass.append(0)
        R_singlepass.append(0)
        N_singlepass.append(0)
        D_singlepass.append(0)
        B_singlepass.append(0)
        C_singlepass.append(0)
        E_singlepass.append(0)
        Q_singlepass.append(0)
        Z_singlepass.append(0)
        G_singlepass.append(0)
        H_singlepass.append(0)
        I_singlepass.append(0)
        L_singlepass.append(0)
        K_singlepass.append(0)
        M_singlepass.append(0)
        F_singlepass.append(0)
        P_singlepass.append(0)
        S_singlepass.append(0)
        T_singlepass.append(0)
        W_singlepass.append(0)
        Y_singlepass.append(0)
        V_singlepass.append(0)

    list_of_tmh_segments = []
    list_of_singlepass_tmh_segments = []

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
            N_flank_sequence = entry[7]
            C_flank_sequence = entry[8]
            tmh_number = entry[9]
            total_tmd_count = entry[10]

            correction_number = 0

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


            if int(total_tmd_count) == 1:
                list_of_singlepass_tmh_segments.append(tmh_segment)




    for i in list_of_singlepass_tmh_segments:

        for position, residue in enumerate(str(i)):
            if residue == "A":
                try:
                    A_singlepass[position] = A_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "R":
                try:
                    R_singlepass[position] = R_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "N":
                try:
                    N_singlepass[position] = N_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "D":
                try:
                    D_singlepass[position] = D_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "B":
                B_singlepass[position] = B_singlepass[position] + 1
            if residue == "C":
                try:
                    C_singlepass[position] = C_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "E":
                try:
                    E_singlepass[position] = E_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "Q":
                try:
                    Q_singlepass[position] = Q_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "Z":
                Z_singlepass[position] = Z_singlepass[position] + 1
            if residue == "G":
                try:
                    G_singlepass[position] = G_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "H":
                try:
                    H_singlepass[position] = H_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "I":
                try:
                    I_singlepass[position] = I_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "L":
                try:
                    L_singlepass[position] = L_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "K":
                try:
                    K_singlepass[position] = K_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "M":
                try:
                    M_singlepass[position] = M_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "F":
                try:
                    F_singlepass[position] = F_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "P":
                try:
                    P_singlepass[position] = P_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "S":
                try:
                    S_singlepass[position] = S_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "T":
                try:
                    T_singlepass[position] = T_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "W":
                try:
                    W_singlepass[position] = W_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "Y":
                try:
                    Y_singlepass[position] = Y_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "V":
                try:
                    V_singlepass[position] = V_singlepass[position] + 1
                except IndexError:
                    pass
            if residue == "J":
                pass


    list(itertools.chain.from_iterable(list_of_singlepass_tmh_segments))



    L_singlepass_normalised = []
    for i in L_singlepass:
        normalised_value = (100 * i)/ sum(L_singlepass)
        L_singlepass_normalised.append(normalised_value)

    R_singlepass_normalised=[]
    for i in R_singlepass:
        normalised_value = (100 * i)/ sum(R_singlepass)
        R_singlepass_normalised.append(normalised_value)

    D_singlepass_normalised=[]
    for i in D_singlepass:
        normalised_value = (100 * i)/ sum(D_singlepass)
        D_singlepass_normalised.append(normalised_value)

    K_singlepass_normalised=[]
    for i in K_singlepass:
        normalised_value = (100 * i)/ sum(K_singlepass)
        K_singlepass_normalised.append(normalised_value)

    E_singlepass_normalised=[]
    for i in E_singlepass:
        normalised_value = (100 * i)/ sum(E_singlepass)
        E_singlepass_normalised.append(normalised_value)

    sequence_position=[]
    for position_number, item in enumerate(L_singlepass):
        sequence_number=position_number -(max_sequence_length/2)
        sequence_position.append(sequence_number)

    values_for_noise_L_inside = L_singlepass_normalised[int((max_sequence_length/2)-30):int((max_sequence_length/2)-25)]
    values_for_noise_L_outside = L_singlepass_normalised[int((max_sequence_length/2)+25):int((max_sequence_length/2)+30)]
    noise_L_inside = np.mean(values_for_noise_L_inside)
    noise_L_outside = np.mean(values_for_noise_L_outside)
    noise_normalised_value_singlepass_L = np.mean([noise_L_inside, noise_L_outside])

    values_for_noise_K_inside = K_singlepass_normalised[int((max_sequence_length/2)-30):int((max_sequence_length/2)-25)]
    values_for_noise_K_outside = K_singlepass_normalised[int((max_sequence_length/2)+25):int((max_sequence_length/2)+30)]
    noise_K_inside = np.mean(values_for_noise_K_inside)
    noise_K_outside = np.mean(values_for_noise_K_outside)
    noise_normalised_value_singlepass_K = np.mean([noise_K_inside, noise_K_outside])

    values_for_noise_R_inside = R_singlepass_normalised[int((max_sequence_length/2)-30):int((max_sequence_length/2)-25)]
    values_for_noise_R_outside = R_singlepass_normalised[int((max_sequence_length/2)+25):int((max_sequence_length/2)+30)]
    noise_R_inside = np.mean(values_for_noise_R_inside)
    noise_R_outside = np.mean(values_for_noise_R_outside)
    noise_normalised_value_singlepass_R = np.mean([noise_R_inside, noise_R_outside])

    values_for_noise_D_inside = D_singlepass_normalised[int((max_sequence_length/2)-30):int((max_sequence_length/2)-25)]
    values_for_noise_D_outside = D_singlepass_normalised[int((max_sequence_length/2)+25):int((max_sequence_length/2)+30)]
    noise_D_inside = np.mean(values_for_noise_D_inside)
    noise_D_outside = np.mean(values_for_noise_D_outside)
    noise_normalised_value_singlepass_D = np.mean([noise_D_inside, noise_D_outside])

    values_for_noise_E_inside = E_singlepass_normalised[int((max_sequence_length/2)-30):int((max_sequence_length/2)-25)]
    values_for_noise_E_outside = E_singlepass_normalised[int((max_sequence_length/2)+25):int((max_sequence_length/2)+30)]
    noise_E_inside = np.mean(values_for_noise_E_inside)
    noise_E_outside = np.mean(values_for_noise_E_outside)
    noise_normalised_value_singlepass_E = np.mean([noise_E_inside, noise_E_outside])


    values_for_flank_E_inside = E_singlepass_normalised[int((max_sequence_length/2)-20):int((max_sequence_length/2)-10)]
    values_for_flank_E_outside = E_singlepass_normalised[int((max_sequence_length/2)+10):int((max_sequence_length/2)+20)]

    noise_normalised_E_inside_singlepass = np.mean(values_for_flank_E_inside)
    noise_normalised_E_outside_singlepass = np.mean(values_for_flank_E_outside)

    noise_normalised_D_inside_singlepass = np.mean(D_singlepass_normalised[int((max_sequence_length/2)-20):int((max_sequence_length/2)-10)])
    noise_normalised_D_outside_singlepass = np.mean(D_singlepass_normalised[int((max_sequence_length/2)+10):int((max_sequence_length/2)+20)])

    noise_normalised_R_inside_singlepass = np.mean(R_singlepass_normalised[int((max_sequence_length/2)-20):int((max_sequence_length/2)-10)])
    noise_normalised_R_outside_singlepass = np.mean(R_singlepass_normalised[int((max_sequence_length/2)+10):int((max_sequence_length/2)+20)])

    noise_normalised_K_inside_singlepass = np.mean(K_singlepass_normalised[int((max_sequence_length/2)-20):int((max_sequence_length/2)-10)])
    noise_normalised_K_outside_singlepass = np.mean(K_singlepass_normalised[int((max_sequence_length/2)+10):int((max_sequence_length/2)+20)])

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

    plt.plot((-50, +50), (noise_normalised_value_singlepass_L, noise_normalised_value_singlepass_L), linestyle='--', color='blue')
    plt.plot((-50, +50), (noise_normalised_value_singlepass_K, noise_normalised_value_singlepass_K), linestyle='--', color='peachpuff')
    plt.plot((-50, +50), (noise_normalised_value_singlepass_R, noise_normalised_value_singlepass_R), linestyle='--', color='orange')
    plt.plot((-50, +50), (noise_normalised_value_singlepass_D, noise_normalised_value_singlepass_D), linestyle='--', color='purple')
    plt.plot((-50, +50), (noise_normalised_value_singlepass_E, noise_normalised_value_singlepass_E), linestyle='--', color='thistle')
    print "noise_normalised_value_singlepass_L", noise_normalised_value_singlepass_L
    print "noise_normalised_value_singlepass_K", noise_normalised_value_singlepass_K
    print "noise_normalised_value_singlepass_R", noise_normalised_value_singlepass_R
    print "noise_normalised_value_singlepass_D", noise_normalised_value_singlepass_D
    print "noise_normalised_value_singlepass_E", noise_normalised_value_singlepass_E

    plt.plot((-20, -10), (noise_normalised_E_inside_singlepass, noise_normalised_E_inside_singlepass), linestyle='-', linewidth=4,  color='thistle')
    plt.plot((10, 20), (noise_normalised_E_outside_singlepass, noise_normalised_E_outside_singlepass), linestyle='-', linewidth=4,  color='thistle')
    print "noise_normalised_E_inside_singlepass", noise_normalised_E_inside_singlepass
    print "noise_normalised_E_outside_singlepass", noise_normalised_E_outside_singlepass

    plt.plot((-20, -10), (noise_normalised_D_inside_singlepass, noise_normalised_D_inside_singlepass), linestyle='-', linewidth=4,  color='purple')
    plt.plot((10, 20), (noise_normalised_D_outside_singlepass, noise_normalised_D_outside_singlepass), linestyle='-', linewidth=4,  color='purple')
    print "noise_normalised_D_inside_singlepass", noise_normalised_D_inside_singlepass
    print "noise_normalised_D_outside_singlepass", noise_normalised_D_outside_singlepass

    plt.plot((-20, -10), (noise_normalised_R_inside_singlepass, noise_normalised_R_inside_singlepass), linestyle='-', linewidth=4,  color='orange')
    plt.plot((10, 20), (noise_normalised_R_outside_singlepass, noise_normalised_R_outside_singlepass), linestyle='-', linewidth=4,  color='orange')
    print "noise_normalised_R_inside_singlepass", noise_normalised_R_inside_singlepass
    print "noise_normalised_R_outside_singlepass", noise_normalised_R_outside_singlepass

    plt.plot((-20, -10), (noise_normalised_K_inside_singlepass, noise_normalised_K_inside_singlepass), linestyle='-', linewidth=4,  color='peachpuff')
    plt.plot((10, 20), (noise_normalised_K_outside_singlepass, noise_normalised_K_outside_singlepass), linestyle='-', linewidth=4,  color='peachpuff')
    print "noise_normalised_K_inside_singlepass", noise_normalised_K_inside_singlepass
    print "noise_normalised_K_outside_singlepass", noise_normalised_K_outside_singlepass


    plt.plot(sequence_position, L_singlepass_normalised, linestyle='-', marker='',  linewidth=2,  color='blue')
    plt.plot(sequence_position, R_singlepass_normalised, linestyle='-', marker='',  linewidth=2,  color='orange')
    plt.plot(sequence_position, K_singlepass_normalised, linestyle='-', marker='',  linewidth=2,  color='peachpuff')
    plt.plot(sequence_position, D_singlepass_normalised, linestyle='-', marker='',  linewidth=2,  color='purple')
    plt.plot(sequence_position, E_singlepass_normalised, linestyle='-', marker='',  linewidth=2,  color='thistle')

    font = {'size':18, 'fontname':'Helvetica Neue Light'}

    plt.xlabel('Sequence Position', **font)
    plt.ylabel('Relative Percentage', **font)
    plt.tick_params(labelsize=16)
    pylab.xlim([-30,30])
    pylab.ylim([0, 6.6])
    timestamp = str(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    filename = timestamp + file.replace('.csv', '') + "_singlepass.pdf"
    plt.gcf().subplots_adjust(bottom=0.2)
    plt.savefig(filename)
    #plt.show()

    plt.clf()
    plt.cla()
    print file
    print "Single-pass:", len(list_of_singlepass_tmh_segments)


for file in list_of_files:
    print "\n",file
    results=[]
    with open(file) as inputfile:
        for line in inputfile:
            results.append(line.strip().split(','))

    total_multipass_inner_insidefraction_sequences = []
    total_multipass_outer_insidefraction_sequences = []

    total_multipass_inner_outsidefraction_sequences = []
    total_multipass_outer_outsidefraction_sequences = []

    list_of_slices=[]

    multipass_helices = 0

    A_multipass = [0]
    R_multipass = [0]
    N_multipass = [0]
    D_multipass = [0]
    B_multipass = [0]
    C_multipass = [0]
    E_multipass = [0]
    Q_multipass = [0]
    Z_multipass = [0]
    G_multipass = [0]
    H_multipass = [0]
    I_multipass = [0]
    L_multipass = [0]
    K_multipass = [0]
    M_multipass = [0]
    F_multipass = [0]
    P_multipass = [0]
    S_multipass = [0]
    T_multipass = [0]
    W_multipass = [0]
    Y_multipass = [0]
    V_multipass = [0]

    maximum_tmd_length = 0
    maximum_segment_length =0
    maxCflank_length=0
    maxNflank_length=0
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
            N_flank_sequence = entry[7]
            C_flank_sequence = entry[8]
            tmh_number = entry[9]
            total_tmd_count = entry[10]

            if "Outside" in str(n_terminal_start):
                inner_flank_length = len(C_flank_sequence)
            if "Inside" in str(n_terminal_start):
                inner_flank_length = len(N_flank_sequence)
            tmd_length=len(tmh_sequence)

            if tmd_length > maximum_tmd_length_for_alignment:
                maximum_tmd_length_for_alignment = tmd_length
            if inner_flank_length > maximum_inner_flank_length:
                maximum_inner_flank_length = inner_flank_length

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

    for _ in range(max_sequence_length):
        A_multipass.append(0)
        R_multipass.append(0)
        N_multipass.append(0)
        D_multipass.append(0)
        B_multipass.append(0)
        C_multipass.append(0)
        E_multipass.append(0)
        Q_multipass.append(0)
        Z_multipass.append(0)
        G_multipass.append(0)
        H_multipass.append(0)
        I_multipass.append(0)
        L_multipass.append(0)
        K_multipass.append(0)
        M_multipass.append(0)
        F_multipass.append(0)
        P_multipass.append(0)
        S_multipass.append(0)
        T_multipass.append(0)
        W_multipass.append(0)
        Y_multipass.append(0)
        V_multipass.append(0)

    list_of_tmh_segments = []
    list_of_multipass_tmh_segments = []
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
            N_flank_sequence = entry[7]
            C_flank_sequence = entry[8]
            tmh_number = entry[9]
            total_tmd_count = entry[10]

            correction_number = 0

            if "Outside" in str(n_terminal_start):
                tmh_unaltered_sequence = str(N_flank_sequence) + str(tmh_sequence) + str(C_flank_sequence)
                tmh_reversed_sequence = tmh_unaltered_sequence[::-1]
                correction_number = ((maximum_tmd_length_for_alignment/2) + maximum_inner_flank_length) - ((len(tmh_sequence)/ 2) + len(C_flank_sequence))
                tmh_segment = "J" * (int(correction_number) + 1)
                tmh_segment = tmh_segment + tmh_reversed_sequence
            if "Inside" in str(n_terminal_start):
                tmh_unaltered_sequence = str(N_flank_sequence) + str(tmh_sequence) + str(C_flank_sequence)
                correction_number = ((maximum_tmd_length_for_alignment/2) + maximum_inner_flank_length) - ((len(tmh_sequence)/ 2) + len(N_flank_sequence))
                tmh_segment = "J" * (int(correction_number) + 1)
                tmh_segment = tmh_segment + tmh_unaltered_sequence
            list_of_tmh_segments.append(tmh_segment)


            if int(total_tmd_count) > 1:
                list_of_multipass_tmh_segments.append(tmh_segment)



    for i in list_of_multipass_tmh_segments:

        for position, residue in enumerate(str(i)):
            if residue == "A":
                try:
                    A_multipass[position] = A_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "R":
                try:
                    R_multipass[position] = R_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "N":
                try:
                    N_multipass[position] = N_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "D":
                try:
                    D_multipass[position] = D_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "B":
                B_multipass[position] = B_multipass[position] + 1
            if residue == "C":
                try:
                    C_multipass[position] = C_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "E":
                try:
                    E_multipass[position] = E_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "Q":
                try:
                    Q_multipass[position] = Q_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "Z":
                Z_multipass[position] = Z_multipass[position] + 1
            if residue == "G":
                try:
                    G_multipass[position] = G_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "H":
                try:
                    H_multipass[position] = H_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "I":
                try:
                    I_multipass[position] = I_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "L":
                try:
                    L_multipass[position] = L_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "K":
                try:
                    K_multipass[position] = K_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "M":
                try:
                    M_multipass[position] = M_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "F":
                try:
                    F_multipass[position] = F_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "P":
                try:
                    P_multipass[position] = P_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "S":
                try:
                    S_multipass[position] = S_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "T":
                try:
                    T_multipass[position] = T_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "W":
                try:
                    W_multipass[position] = W_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "Y":
                try:
                    Y_multipass[position] = Y_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "V":
                try:
                    V_multipass[position] = V_multipass[position] + 1
                except IndexError:
                    pass
            if residue == "J":
                pass


    list(itertools.chain.from_iterable(list_of_multipass_tmh_segments))

    L_multipass_normalised = []
    for i in L_multipass:
        normalised_value = (100 * i)/ sum(L_multipass)
        L_multipass_normalised.append(normalised_value)

    R_multipass_normalised=[]
    for i in R_multipass:
        normalised_value = (100 * i)/ sum(R_multipass)
        R_multipass_normalised.append(normalised_value)

    D_multipass_normalised=[]
    for i in D_multipass:
        normalised_value = (100 * i)/ sum(D_multipass)
        D_multipass_normalised.append(normalised_value)

    K_multipass_normalised=[]
    for i in K_multipass:
        normalised_value = (100 * i)/ sum(K_multipass)
        K_multipass_normalised.append(normalised_value)

    E_multipass_normalised=[]
    for i in E_multipass:
        normalised_value = (100 * i)/ sum(E_multipass)
        E_multipass_normalised.append(normalised_value)

    sequence_position=[]
    for position_number, item in enumerate(L_multipass):
        sequence_number=position_number -(max_sequence_length/2)
        sequence_position.append(sequence_number)

    values_for_noise_L_inside = L_multipass_normalised[int((max_sequence_length/2)-30):int((max_sequence_length/2)-25)]
    values_for_noise_L_outside = L_multipass_normalised[int((max_sequence_length/2)+25):int((max_sequence_length/2)+30)]
    noise_L_inside = np.mean(values_for_noise_L_inside)
    noise_L_outside = np.mean(values_for_noise_L_outside)
    noise_normalised_value_multipass_L = np.mean([noise_L_inside, noise_L_outside])

    values_for_noise_K_inside = K_multipass_normalised[int((max_sequence_length/2)-30):int((max_sequence_length/2)-25)]
    values_for_noise_K_outside = K_multipass_normalised[int((max_sequence_length/2)+25):int((max_sequence_length/2)+30)]
    noise_K_inside = np.mean(values_for_noise_K_inside)
    noise_K_outside = np.mean(values_for_noise_K_outside)
    noise_normalised_value_multipass_K = np.mean([noise_K_inside, noise_K_outside])

    values_for_noise_R_inside = R_multipass_normalised[int((max_sequence_length/2)-30):int((max_sequence_length/2)-25)]
    values_for_noise_R_outside = R_multipass_normalised[int((max_sequence_length/2)+25):int((max_sequence_length/2)+30)]
    noise_R_inside = np.mean(values_for_noise_R_inside)
    noise_R_outside = np.mean(values_for_noise_R_outside)
    noise_normalised_value_multipass_R = np.mean([noise_R_inside, noise_R_outside])

    values_for_noise_D_inside = D_multipass_normalised[int((max_sequence_length/2)-30):int((max_sequence_length/2)-25)]
    values_for_noise_D_outside = D_multipass_normalised[int((max_sequence_length/2)+25):int((max_sequence_length/2)+30)]
    noise_D_inside = np.mean(values_for_noise_D_inside)
    noise_D_outside = np.mean(values_for_noise_D_outside)
    noise_normalised_value_multipass_D = np.mean([noise_D_inside, noise_D_outside])

    values_for_noise_E_inside = E_multipass_normalised[int((max_sequence_length/2)-30):int((max_sequence_length/2)-25)]
    values_for_noise_E_outside = E_multipass_normalised[int((max_sequence_length/2)+25):int((max_sequence_length/2)+30)]
    noise_E_inside = np.mean(values_for_noise_E_inside)
    noise_E_outside = np.mean(values_for_noise_E_outside)
    noise_normalised_value_multipass_E = np.mean([noise_E_inside, noise_E_outside])


    values_for_flank_E_inside = E_multipass_normalised[int((max_sequence_length/2)-20):int((max_sequence_length/2)-10)]
    values_for_flank_E_outside = E_multipass_normalised[int((max_sequence_length/2)+10):int((max_sequence_length/2)+20)]

    noise_normalised_E_inside_multipass = np.mean(values_for_flank_E_inside)
    noise_normalised_E_outside_multipass = np.mean(values_for_flank_E_outside)

    noise_normalised_D_inside_multipass = np.mean(D_multipass_normalised[int((max_sequence_length/2)-20):int((max_sequence_length/2)-10)])
    noise_normalised_D_outside_multipass = np.mean(D_multipass_normalised[int((max_sequence_length/2)+10):int((max_sequence_length/2)+20)])

    noise_normalised_R_inside_multipass = np.mean(R_multipass_normalised[int((max_sequence_length/2)-20):int((max_sequence_length/2)-10)])
    noise_normalised_R_outside_multipass = np.mean(R_multipass_normalised[int((max_sequence_length/2)+10):int((max_sequence_length/2)+20)])

    noise_normalised_K_inside_multipass = np.mean(K_multipass_normalised[int((max_sequence_length/2)-20):int((max_sequence_length/2)-10)])
    noise_normalised_K_outside_multipass = np.mean(K_multipass_normalised[int((max_sequence_length/2)+10):int((max_sequence_length/2)+20)])

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

    plt.plot((-50, +50), (noise_normalised_value_multipass_L, noise_normalised_value_multipass_L), linestyle='--', color='blue')
    plt.plot((-50, +50), (noise_normalised_value_multipass_K, noise_normalised_value_multipass_K), linestyle='--', color='peachpuff')
    plt.plot((-50, +50), (noise_normalised_value_multipass_R, noise_normalised_value_multipass_R), linestyle='--', color='orange')
    plt.plot((-50, +50), (noise_normalised_value_multipass_D, noise_normalised_value_multipass_D), linestyle='--', color='purple')
    plt.plot((-50, +50), (noise_normalised_value_multipass_E, noise_normalised_value_multipass_E), linestyle='--', color='thistle')

    print "noise_normalised_value_multipass_L", noise_normalised_value_multipass_L
    print "noise_normalised_value_multipass_K", noise_normalised_value_multipass_K
    print "noise_normalised_value_multipass_R", noise_normalised_value_multipass_R
    print "noise_normalised_value_multipass_D", noise_normalised_value_multipass_D
    print "noise_normalised_value_multipass_E", noise_normalised_value_multipass_E

    plt.plot((-20, -10), (noise_normalised_E_inside_multipass, noise_normalised_E_inside_multipass), linestyle='-', linewidth=4,  color='thistle')
    plt.plot((10, 20), (noise_normalised_E_outside_multipass, noise_normalised_E_outside_multipass), linestyle='-', linewidth=4,  color='thistle')
    print "noise_normalised_E_inside_multipass", noise_normalised_E_inside_multipass
    print "noise_normalised_E_outside_multipass", noise_normalised_E_outside_multipass

    plt.plot((-20, -10), (noise_normalised_D_inside_multipass, noise_normalised_D_inside_multipass), linestyle='-', linewidth=4,  color='purple')
    plt.plot((10, 20), (noise_normalised_D_outside_multipass, noise_normalised_D_outside_multipass), linestyle='-', linewidth=4,  color='purple')
    print "noise_normalised_D_inside_multipass", noise_normalised_D_inside_multipass
    print "noise_normalised_D_outside_multipass", noise_normalised_D_outside_multipass

    plt.plot((-20, -10), (noise_normalised_R_inside_multipass, noise_normalised_R_inside_multipass), linestyle='-', linewidth=4,  color='orange')
    plt.plot((10, 20), (noise_normalised_R_outside_multipass, noise_normalised_R_outside_multipass), linestyle='-', linewidth=4,  color='orange')
    print "noise_normalised_R_inside_multipass", noise_normalised_R_inside_multipass
    print "noise_normalised_R_outside_multipass", noise_normalised_R_outside_multipass

    plt.plot((-20, -10), (noise_normalised_K_inside_multipass, noise_normalised_K_inside_multipass), linestyle='-', linewidth=4,  color='peachpuff')
    plt.plot((10, 20), (noise_normalised_K_outside_multipass, noise_normalised_K_outside_multipass), linestyle='-', linewidth=4,  color='peachpuff')
    print "noise_normalised_K_inside_multipass", noise_normalised_K_inside_multipass
    print "noise_normalised_K_outside_multipass", noise_normalised_K_outside_multipass


    plt.plot(sequence_position, L_multipass_normalised, linestyle='-', marker='',  linewidth=2,  color='blue')
    plt.plot(sequence_position, R_multipass_normalised, linestyle='-', marker='',  linewidth=2,  color='orange')
    plt.plot(sequence_position, K_multipass_normalised, linestyle='-', marker='',  linewidth=2,  color='peachpuff')
    plt.plot(sequence_position, D_multipass_normalised, linestyle='-', marker='',  linewidth=2,  color='purple')
    plt.plot(sequence_position, E_multipass_normalised, linestyle='-', marker='',  linewidth=2,  color='thistle')

    font = {'size':18, 'fontname':'Helvetica Neue Light'}

    plt.xlabel('Sequence Position', **font)
    plt.ylabel('Relative Percentage', **font)
    plt.tick_params(labelsize=16)
    pylab.xlim([-30,30])
    pylab.ylim([0, 6])
    timestamp = str(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    filename = timestamp + file.replace('.csv', '') + "_multipass.pdf"
    plt.gcf().subplots_adjust(bottom=0.2)
    plt.savefig(filename)
    #plt.show()

    plt.clf()
    plt.cla()
    print file
    print "Multi-pass:", len(list_of_multipass_tmh_segments)
