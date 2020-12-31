"""
This file is designed for a project to integrate
the GOR protein secondary structure prediction algorithm 
with the TALOS algorithm. The approx_gor3, approx_gor4, approx_gor5
functions implement GOR 3, 4 and 5 algorithms respectively

author: Parsia Basimfar.
"""

import math
import pandas as pd
import numpy as np

pd.options.mode.chained_assignment = None
df = pd.read_excel('GOR.xlsx', convert_float=True)
df = df.replace(np.nan, 'NaN', regex=True)
residue = df['RESIDUE']


def print_full(x):
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')


def convert_col(col_name):  # converts floats to ints in a given column
    col = df[col_name]
    for i in range(len(col)):
        if isinstance(col[i], float) is True:
            col[i] = int(col[i])
        else:
            continue
    return df


def convert_dsspnames(dataframe):  # converts the 8-letter system to 3-letter system
    col = dataframe['STRUCTURE']
    for i in range(len(col)):
        if col[i] == 'G' or col[i] == 'I':
            col[i] = 'H'
        elif col[i] == 'B':
            col[i] = 'E'
        elif col[i] == 'T' or col[i] == 'S':
            col[i] = 'C'
        else:
            pass


# count NaN residues in database
NaN = []
for j in range(len(df)):
    if 'NaN' in df.iloc[j, 3]:
        NaN.append(1)

amino_acids = {'A': 'alanin', 'R': 'arginine', 'N': 'asparagine', 'D': 'aspartic acid',
               'C': 'cystein', 'E': 'glutamic acid', 'Q': 'glutamine', 'G': 'glycine',
               'H': 'histidine', 'I': 'isoleucine', 'L': 'leucine', 'K': 'lysine',
               'M': 'methionine', 'F': 'phenylalanine', 'P': 'proline', 'S': 'serine',
               'T': 'threonine', 'W': 'tryptophan', 'Y': 'tyrosine', 'V': 'valine'}


# Find the Fano Function for each amino acid
def mutual_information(aa):
    convert_dsspnames(df)

    conf_h = []
    conf_e = []
    conf_c = []
    s_h = []
    s_e = []
    s_c = []

    # check this code if it can become more compact
    for i in range(1, len(df)):
        if 'H' in df.iloc[i, 3]:
            s_h.append(1)
        elif 'E' in df.iloc[i, 3]:
            s_e.append(1)
        elif 'C' in df.iloc[i, 3]:
            s_c.append(1)
        else:
            pass

        if aa in df.iloc[i, 2]:
            if df.iloc[i, 3] is not 'NaN':
                if df.iloc[i, 3] == 'H':
                    conf_h.append(1)
                elif df.iloc[i, 3] == 'E':
                    conf_e.append(1)
                else:
                    conf_c.append(1)
    sum_aa = sum(conf_h) + sum(conf_e) + sum(conf_c)

    # Information difference
    fano_h = math.log(sum(conf_h) / (sum(conf_c) + sum(conf_e))) + math.log((sum(s_c) + sum(s_e)) / sum(s_h))
    fano_e = math.log(sum(conf_e) / (sum(conf_c) + sum(conf_h))) + math.log((sum(s_c) + sum(s_h)) / sum(s_e))
    fano_c = math.log(sum(conf_c) / (sum(conf_h) + sum(conf_e))) + math.log((sum(s_h) + sum(s_e)) / sum(s_c))

    if max(fano_c, fano_e, fano_h) == fano_c:
        print('Coil', aa)
    elif max(fano_c, fano_e, fano_h) == fano_e:
        print('Sheet', aa)
    else:
        print('Helix', aa)

    return [fano_h, fano_e, fano_c, aa]


def approx_gor3(target):
    keys = ['A', 'V', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y']
    # Pre-allocation
    mu_info = {}
    for key in keys:
        mu_info[key] = mutual_information(key)

    info_window_h = [[] for _ in range(len(target))]
    info_window_e = [[] for _ in range(len(target))]
    info_window_c = [[] for _ in range(len(target))]
    for i in range(len(target)):
        if 0 <= i < 8:  # ignore the 8 first and 8 last residues(else block)
            info_window_h[i].append(0)
            info_window_e[i].append(0)
            info_window_c[i].append(0)
        elif 8 <= i <= len(target) - 9:
            # here the window is constructed
            for window in range(1, 9):
                info_window_h[i].append(mu_info[target[i]][0] + mu_info[target[i + window]][0])
                info_window_h[i].append(mu_info[target[i - window]][0])
                info_window_e[i].append(mu_info[target[i]][1] + mu_info[target[i + window]][1])
                info_window_e[i].append(mu_info[target[i - window]][1])
                info_window_c[i].append(mu_info[target[i]][2] + mu_info[target[i + window]][2])
                info_window_c[i].append(mu_info[target[i - window]][2])
        else:
            info_window_h[i].append(0)
            info_window_e[i].append(0)
            info_window_c[i].append(0)

    approx = [[] for _ in range(len(target))]
    for i in range(len(target)):
        approx[i].append(sum(info_window_h[i]))
        approx[i].append(sum(info_window_e[i]))
        approx[i].append(sum(info_window_c[i]))

        if approx[i][0] == 0 and approx[i][1] == 0 and approx[i][2] == 0:
            print("%s : -" % target[i])
        elif max(approx[i]) == approx[i][0]:
            print("%s : H" % target[i])
        elif max(approx[i]) == approx[i][1]:
            print("%s : E" % target[i])
        elif max(approx[i]) == approx[i][2]:
            print("%s : C" % target[i])


def approx_gor4(target):
    convert_dsspnames(df)
    keys = ['A', 'V', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y']

    # Pre-allocation
    mu_info = {}
    for key in keys:
        mu_info[key] = mutual_information(key)

    s_h = []
    s_e = []
    s_c = []

    for i in range(1, len(df)):
        if 'H' in df.iloc[i, 3]:
            s_h.append(1)
        elif 'E' in df.iloc[i, 3]:
            s_e.append(1)
        elif 'C' in df.iloc[i, 3]:
            s_c.append(1)
        else:
            pass

    prob_h = sum(s_h) / (sum(s_e) + sum(s_c))  # P(S)/P(n-S)
    prob_e = sum(s_e) / (sum(s_h) + sum(s_c))
    prob_c = sum(s_c) / (sum(s_e) + sum(s_h))

    info_window_h = [[] for _ in range(len(target))]
    info_window_e = [[] for _ in range(len(target))]
    info_window_c = [[] for _ in range(len(target))]

    for i in range(len(target)):
        if 0 <= i < 8:  # ignore the 8 first and 8 last residues(else block)
            info_window_h[i].append(0)
            info_window_e[i].append(0)
            info_window_c[i].append(0)
        elif 8 <= i <= len(target) - 9:
            for window_l in range(-8, 8):
                try:
                    for k in range(7 - window_l):  # 7-j to avoid self-summation

                        info_window_h[i].append((2 / 17) * (math.log(mu_info[target[i + window_l]][0] +
                                                                     mu_info[target[i + window_l + k]][0]) + math.log(
                            prob_h)))

                        info_window_e[i].append((2 / 17) * (math.log(mu_info[target[i + window_l]][1] +
                                                                     mu_info[target[i + window_l + k]][1]) + math.log(
                            prob_e)))

                        info_window_c[i].append((2 / 17) * (math.log(mu_info[target[i + window_l]][2] +
                                                                     mu_info[target[i + window_l + k]][2]) + math.log(
                            prob_c)))
                # ignore 0 probability in the database
                except ValueError:
                    pass
            try:
                for window_r in range(1, 9):
                    info_window_h[i].append(
                        -(15 / 17) * (mu_info[target[i]][0] + mu_info[target[i + window_r]][0] + math.log(
                            prob_h)))
                    info_window_h[i].append(-(15 / 17) * (mu_info[target[i - window_r]][0] + math.log(
                        prob_h)))
                    info_window_e[i].append(
                        -(15 / 17) * (mu_info[target[i]][1] + mu_info[target[i + window_r]][1] + math.log(
                            prob_h)))
                    info_window_e[i].append(-(15 / 17) * (mu_info[target[i - window_r]][1] + math.log(
                        prob_h)))
                    info_window_c[i].append(
                        -(15 / 17) * (mu_info[target[i]][2] + mu_info[target[i + window_r]][2] + math.log(
                            prob_h)))
                    info_window_c[i].append(-(15 / 17) * (mu_info[target[i - window_r]][2] + math.log(
                        prob_h)))
            except ValueError:
                pass
        else:
            info_window_h[i].append(0)
            info_window_e[i].append(0)
            info_window_c[i].append(0)

    approx = [[] for _ in range(len(target))]
    for i in range(len(target)):
        approx[i].append(sum(info_window_h[i]))
        approx[i].append(sum(info_window_e[i]))
        approx[i].append(sum(info_window_c[i]))

        if approx[i][0] == 0 and approx[i][1] == 0 and approx[i][2] == 0:
            print("%s : -" % target[i])
        elif max(approx[i]) == approx[i][0]:
            print("%s : H" % target[i])
        elif max(approx[i]) == approx[i][1]:
            print("%s : E" % target[i])
        elif max(approx[i]) == approx[i][2]:
            print("%s : C" % target[i])


def approx_gor5(target):

    if type(target) != str:
        raise TypeError('The input should be a string of amino acid sequences.')
    if len(target) <= 7:
        raise ValueError('The sequence is too short!')

    print("The length of the query protein is %i amino acids." % len(target))
    convert_dsspnames(df)
    keys = ['A', 'V', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y']

    # Pre-allocation
    mu_info = {}
    for key in keys:
        mu_info[key] = mutual_information(key)

    s_h = []
    s_e = []
    s_c = []

    for i in range(1, len(df)):
        if 'H' in df.iloc[i, 3]:
            s_h.append(1)
        elif 'E' in df.iloc[i, 3]:
            s_e.append(1)
        elif 'C' in df.iloc[i, 3]:
            s_c.append(1)
        else:
            pass

    prob_h = sum(s_h) / (sum(s_e) + sum(s_c))  # P(S)/P(n-S)
    prob_e = sum(s_e) / (sum(s_h) + sum(s_c))
    prob_c = sum(s_c) / (sum(s_e) + sum(s_h))

    info_window_h = [[] for _ in range(len(target))]
    info_window_e = [[] for _ in range(len(target))]
    info_window_c = [[] for _ in range(len(target))]

    if len(target) >= 100:
        d = 6
    elif 51 <= len(target) < 100:
        d = 5
    elif 25 < len(target) < 51:
        d = 4
    elif 7 < len(target) <= 25:
        d = 3

    for i in range(len(target)):
        if 0 <= i < d:  # ignore the 8 first and 8 last residues(else block)
            info_window_h[i].append(0)
            info_window_e[i].append(0)
            info_window_c[i].append(0)
        elif d <= i <= len(target) - d - 1:
            for window_l in range(-d, d):
                try:
                    for k in range(d - 1 - window_l):  # 7-j to avoid self-summation

                        info_window_h[i].append((2 / 2*d + 1) * (math.log(mu_info[target[i + window_l]][0] +
                                                                     mu_info[target[i + window_l + k]][0]) + math.log(
                            prob_h)))

                        info_window_e[i].append((2 / 2*d + 1) * (math.log(mu_info[target[i + window_l]][1] +
                                                                     mu_info[target[i + window_l + k]][1]) + math.log(
                            prob_e)))

                        info_window_c[i].append((2 / 2*d + 1) * (math.log(mu_info[target[i + window_l]][2] +
                                                                     mu_info[target[i + window_l + k]][2]) + math.log(
                            prob_c)))
                # ignore 0 probability in the database
                except ValueError:
                    pass
            try:
                for window_r in range(1, d + 1):
                    info_window_h[i].append(
                        -(2*d - 1 / 2*d + 1) * (mu_info[target[i]][0] + mu_info[target[i + window_r]][0] + math.log(
                            prob_h)))
                    info_window_h[i].append(-(2*d - 1 / 2*d + 1) * (mu_info[target[i - window_r]][0] + math.log(
                        prob_h)))
                    info_window_e[i].append(
                        -(2*d - 1 / 2*d + 1) * (mu_info[target[i]][1] + mu_info[target[i + window_r]][1] + math.log(
                            prob_h)))
                    info_window_e[i].append(-(2*d - 1 / 2*d + 1) * (mu_info[target[i - window_r]][1] + math.log(
                        prob_h)))
                    info_window_c[i].append(
                        -(2*d - 1 / 2*d + 1) * (mu_info[target[i]][2] + mu_info[target[i + window_r]][2] + math.log(
                            prob_h)))
                    info_window_c[i].append(-(2*d - 1 / 2*d + 1) * (mu_info[target[i - window_r]][2] + math.log(
                        prob_h)))
            except ValueError:
                pass
        else:
            info_window_h[i].append(0)
            info_window_e[i].append(0)
            info_window_c[i].append(0)

    approx = [[] for _ in range(len(target))]

    for i in range(len(target)):
        approx[i].append(sum(info_window_h[i]))
        approx[i].append(sum(info_window_e[i]))
        approx[i].append(sum(info_window_c[i]))
        if approx[i][0] == 0 and approx[i][1] == 0 and approx[i][2] == 0:
            print("%s : -" % target[i])
        elif max(approx[i]) == approx[i][0]:
            print("%s : H" % target[i])
        elif max(approx[i]) == approx[i][1]:
            print("%s : E" % target[i])
        elif max(approx[i]) == approx[i][2]:
            print("%s : C" % target[i])


test = "MADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMAR" \
       "KMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK"

