from __future__ import division

list_of_files = [
    "TopDB_5_flanklength_flankclashTrue.csv",
    "Unihuman_5_flanklength_flankclashTrue.csv",
    "TopDB_20_flanklength_flankclashTrue.csv",
    "Unihuman_20_flanklength_flankclashTrue.csv",
    ]

for file in list_of_files:
    A_multipass= 0
    R_multipass= 0
    N_multipass= 0
    D_multipass= 0
    B_multipass= 0
    C_multipass= 0
    E_multipass= 0
    Q_multipass= 0
    Z_multipass= 0
    G_multipass= 0
    H_multipass= 0
    I_multipass= 0
    L_multipass= 0
    K_multipass= 0
    M_multipass= 0
    F_multipass= 0
    P_multipass= 0
    S_multipass= 0
    T_multipass= 0
    W_multipass= 0
    Y_multipass= 0
    V_multipass= 0
    U_multipass=0
    A_singlepass= 0
    R_singlepass= 0
    N_singlepass= 0
    D_singlepass= 0
    B_singlepass= 0
    C_singlepass= 0
    E_singlepass= 0
    Q_singlepass= 0
    Z_singlepass= 0
    G_singlepass= 0
    H_singlepass= 0
    I_singlepass= 0
    L_singlepass= 0
    K_singlepass= 0
    M_singlepass= 0
    F_singlepass= 0
    P_singlepass= 0
    S_singlepass= 0
    T_singlepass= 0
    W_singlepass= 0
    Y_singlepass= 0
    V_singlepass= 0
    U_singlepass= 0

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
            N_flank_sequence = entry[7]
            C_flank_sequence = entry[8]
            tmh_number = entry[9]
            total_tmd_count = entry[10]

            correction_number = 0
            tmh_sequence = str(N_flank_sequence) + str(tmh_sequence) + str(C_flank_sequence)

            if int(total_tmd_count) == 1:
                list_of_singlepass_tmh_segments.append(tmh_sequence)
            elif int(total_tmd_count) > 1:
                list_of_multipass_tmh_segments.append(tmh_sequence)

    for i in list_of_multipass_tmh_segments:

        for position, residue in enumerate(str(i)):

            if residue == "A":

                A_multipass = A_multipass + 1

            if residue == "R":

                R_multipass = R_multipass + 1

            if residue == "N":

                N_multipass = N_multipass + 1

            if residue == "D":

                D_multipass = D_multipass + 1

            if residue == "B":

                B_multipass = B_multipass + 1

            if residue == "C":

                C_multipass = C_multipass + 1

            if residue == "E":

                E_multipass = E_multipass + 1

            if residue == "Q":

                Q_multipass = Q_multipass + 1

            if residue == "Z":

                Z_multipass = Z_multipass + 1

            if residue == "G":

                G_multipass = G_multipass + 1

            if residue == "H":

                H_multipass = H_multipass + 1

            if residue == "I":

                I_multipass = I_multipass + 1

            if residue == "L":

                L_multipass = L_multipass + 1

            if residue == "K":
                K_multipass = K_multipass + 1

            if residue == "M":

                M_multipass = M_multipass + 1

            if residue == "F":

                F_multipass = F_multipass + 1

            if residue == "P":

                P_multipass = P_multipass + 1

            if residue == "S":

                S_multipass = S_multipass + 1

            if residue == "T":

                T_multipass = T_multipass + 1

            if residue == "W":

                W_multipass = W_multipass + 1

            if residue == "Y":

                Y_multipass = Y_multipass + 1

            if residue == "V":

                V_multipass = V_multipass + 1
            if residue == "U":

                U_multipass = U_multipass + 1


            if residue == "J":

                pass

    for i in list_of_singlepass_tmh_segments:

        for position, residue in enumerate(str(i)):

            if residue == "A":

                A_singlepass = A_singlepass + 1

            if residue == "R":

                R_singlepass = R_singlepass + 1

            if residue == "N":

                N_singlepass = N_singlepass + 1

            if residue == "D":

                D_singlepass = D_singlepass + 1

            if residue == "B":

                B_singlepass = B_singlepass + 1

            if residue == "C":

                C_singlepass = C_singlepass + 1

            if residue == "E":

                E_singlepass = E_singlepass + 1

            if residue == "Q":

                Q_singlepass = Q_singlepass + 1

            if residue == "Z":

                Z_singlepass = Z_singlepass + 1

            if residue == "G":

                G_singlepass = G_singlepass + 1

            if residue == "H":

                H_singlepass = H_singlepass + 1

            if residue == "I":

                I_singlepass = I_singlepass + 1

            if residue == "L":

                L_singlepass = L_singlepass + 1

            if residue == "K":
                K_singlepass = K_singlepass + 1

            if residue == "M":

                M_singlepass = M_singlepass + 1

            if residue == "F":

                F_singlepass = F_singlepass + 1

            if residue == "P":

                P_singlepass = P_singlepass + 1

            if residue == "S":

                S_singlepass = S_singlepass + 1

            if residue == "T":

                T_singlepass = T_singlepass + 1

            if residue == "W":

                W_singlepass = W_singlepass + 1

            if residue == "Y":

                Y_singlepass = Y_singlepass + 1

            if residue == "V":

                V_singlepass = V_singlepass + 1

            if residue == "U":

                U_singlepass = U_singlepass + 1

            if residue == "J":

                pass


    print file
    print "Single-pass,", len(list_of_singlepass_tmh_segments)
    print "A,", A_singlepass
    print "R,", R_singlepass
    print "N,", N_singlepass
    print "D,", D_singlepass
    print "B,", B_singlepass
    print "C,", C_singlepass
    print "E,", E_singlepass
    print "Q,", Q_singlepass
    print "Z,", Z_singlepass
    print "G,", G_singlepass
    print "H,", H_singlepass
    print "I,", I_singlepass
    print "L,", L_singlepass
    print "K,", K_singlepass
    print "M,", M_singlepass
    print "F,", F_singlepass
    print "P,", P_singlepass
    print "S,", S_singlepass
    print "T,", T_singlepass
    print "W,", W_singlepass
    print "Y,", Y_singlepass
    print "V,", V_singlepass
    print "U,", U_singlepass


    print "Multi-pass,", len(list_of_multipass_tmh_segments)
    print "A,", A_multipass
    print "R,", R_multipass
    print "N,", N_multipass
    print "D,", D_multipass
    print "B,", B_multipass
    print "C,", C_multipass
    print "E,", E_multipass
    print "Q,", Q_multipass
    print "Z,", Z_multipass
    print "G,", G_multipass
    print "H,", H_multipass
    print "I,", I_multipass
    print "L,", L_multipass
    print "K,", K_multipass
    print "M,", M_multipass
    print "F,", F_multipass
    print "P,", P_multipass
    print "S,", S_multipass
    print "T,", T_multipass
    print "W,", W_multipass
    print "Y,", Y_multipass
    print "V,", V_multipass
    print "U,", U_multipass
