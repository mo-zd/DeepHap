import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy
import time
from setting import BAF_ASSIGNMENT_DICT
from setting import LOGR_ASSIGNMENT_DICT
from setting import ANOMALY_DICT
from plot import plot_baf
from plot import plot_logr
from read_write import read_baf_file
from read_write import read_logr_file
import pprint
# import ruptures as rpt
import os.path as osp
from collections import Counter
import random


def denoise(array, allowed_values, window_size=10):
    starts = [0] * int(window_size / 2) + list(range(0, len(array) - int(window_size / 2)))
    ends = list(range(int(window_size / 2), len(array))) + [len(array) - 1] * int(window_size / 2)
    windows = list(zip(starts, ends))
    out = []
    for window in windows:
        counts = Counter(array[window[0]: window[1]])
        val = allowed_values[np.argmax([counts[i] for i in allowed_values])]
        out.append(val)
    return out


def process_baf_data(baf_data, window_size=0):
    global BAF_ASSIGNMENT_DICT
    baf_info = {}
    for ind, (chr_name, loc, baf) in baf_data.items():
        baf_info[ind] = (chr_name, loc, BAF_ASSIGNMENT_DICT[int(baf * 12.)])
    if window_size:
        out = sorted(list(baf_info.items()), key=lambda x: x[1][1])
        keys, vals = list(zip(*out))
        chr_names, locs, bafs = zip(*vals)
        allowed_values = sorted(tuple(set(bafs)))
        bafs = denoise(bafs, allowed_values, window_size)
        baf_info = dict(zip(keys, list(zip(chr_names, locs, bafs))))
    return baf_info


def process_logr_data(logr_data, window_size=0):
    global LOGR_ASSIGNMENT_DICT
    logr_info = {}
    for ind, (chr_name, loc, logr) in logr_data.items():
        logr_info[ind] = (chr_name, loc, LOGR_ASSIGNMENT_DICT[int((logr + 2.) * 2.)])
    if window_size:
        out = sorted(list(logr_info.items()), key=lambda x: x[1][1])
        keys, vals = list(zip(*out))
        chr_names, locs, logrs = zip(*vals)
        allowed_values = sorted(tuple(set(logrs)))
        logrs = denoise(logrs, allowed_values, window_size)
        logr_info = dict(zip(keys, list(zip(chr_names, locs, logrs))))

    return logr_info


def nearest1(query, pivot):
    out = []
    lb, ub = -1000, query[0]
    mid = (lb + ub)/2.
    query = query[1:]
    for el in pivot:
        while el >= ub:
            lb = ub
            if query:
                ub = query.pop(0)
            else:
                break
        mid = (lb + ub)/2.
        if el < mid:
            out.append(lb)
        else:
            out.append(ub)
    return out


def nearest(query, pivot):
    out = []
    lb, ub = -1000, query[0][1]
    indl, indu = query[0][0], query[0][0]
    mid = (lb + ub)/2.
    query = query[1:]
    for el in pivot:
        while el[1] >= ub:
            lb = ub
            indl = indu
            if query:
                indu, ub = query.pop(0)
            else:
                break
        mid = (lb + ub)/2.
        if el[1] < mid:
            out.append(indl)
        else:
            out.append(indu)
    return out


def generate_neighbor_locs(m1_info, m2_info, p1_info, p2_info, logr_info):
    m1 = [(k, v[1]) for k, v in m1_info.items()]
    m2 = [(k, v[1]) for k, v in m2_info.items()]
    p1 = [(k, v[1]) for k, v in p1_info.items()]
    p2 = [(k, v[1]) for k, v in p2_info.items()]
    logr = [(k, v[1]) for k, v in logr_info.items()]

    m1_list = [el[0] for el in m1]
    m2_list = [el[0] for el in m2]
    p1_list = [el[0] for el in p1]
    p2_list = [el[0] for el in p2]
    logr_list = [el[0] for el in logr]

    set1 = set(zip(m1_list, nearest(m2, m1), nearest(p1, m1), nearest(p2, m1), nearest(logr, m1)))
    set2 = set(zip(nearest(m1, m2), m2_list, nearest(p1, m2), nearest(p2, m2), nearest(logr, m2)))
    set3 = set(zip(nearest(m1, p1), nearest(m2, p1), p1_list, nearest(p2, p1), nearest(logr, p1)))
    set4 = set(zip(nearest(m1, p2), nearest(m2, p2), nearest(p1, p2), p2_list, nearest(logr, p2)))
    set5 = set(zip(nearest(m1, logr), nearest(m2, logr), nearest(p1, logr), nearest(p2, logr), logr_list))
    out = set1.union(set2, set3, set4, set5)
    return out


def generate_hap_signatures(m1_info, m2_info, p1_info, p2_info, logr_info):
    out = generate_neighbor_locs(m1_info, m2_info, p1_info, p2_info, logr_info)
    out = {(e1, e2, e3, e4, e5): (m1_info[e1][2], m2_info[e2][2], p1_info[e3][2], p2_info[e4][2], logr_info[e5][2]) for
     e1, e2, e3, e4, e5 in out}
    return out


def assign_anomaly(hap_signatures):
    return {k: ANOMALY_DICT[v] for k, v in hap_signatures.items() if v in ANOMALY_DICT.keys()}


def generate_genome_logr(file):
    df = pd.read_csv(file, sep='\t')
    chr_list = set(df['Chr'].values)
    if 'X' in chr_list:
        chr_list = [str(i) for i in range(1, 23)].append('X')
    logr_info_list = []
    for chr_name in chr_list:
        df1 = df.loc[df['Chr']==chr_name]
        logr_info_list.append((chr_name, process_logr_data(df1)))


if __name__ == '__main__':
     # start = time.time()
     # logr_data = read_logr_file('./files/E03_Bl004_LogRsAvgWindow_1.txt')
     # logr_info = process_logr_data(logr_data)
     # plot_logr(logr_data, logr_info)
     # print('timeing is:', time.time()-start)
     # start = time.time()
     # baf_data = read_baf_file('./files/E03_Bl004_M1_Seg_1.txt')
     # baf_info = process_baf_data(baf_data)
     # plot_baf(baf_data, baf_info)
     # print('timeing is:', time.time()-start)

    #=====================================================

     start = time.time()
     window_size = 600
     m1_data = read_baf_file('./files/E06_Bl001/E08_Bl001_M1.txt')
     m1_info = process_baf_data(m1_data, window_size=window_size)
     m2_data = read_baf_file('./files/E06_Bl001/E08_Bl001_M2.txt')
     m2_info = process_baf_data(m2_data, window_size=window_size)
     p1_data = read_baf_file('./files/E06_Bl001/E08_Bl001_P1.txt')
     p1_info = process_baf_data(p1_data, window_size=window_size)
     p2_data = read_baf_file('./files/E06_Bl001/E08_Bl001_P2.txt')
     p2_info = process_baf_data(p2_data, window_size=window_size)
     logr_data = read_logr_file('./files/E06_Bl001/E08_Bl001_LOGR.txt')
     logr_info = process_logr_data(logr_data, window_size=window_size)

     plot_logr(logr_data, logr_info)
     plot_baf(m1_data, m1_info)
     plot_baf(m2_data, m2_info)
     plot_baf(p1_data, p1_info)
     plot_baf(p2_data, p2_info)


     out = generate_hap_signatures(m1_info, m2_info, p1_info, p2_info, logr_info)
     print('timing:', time.time() - start)
     print('number of elements of out is', len(out))
     anomaly_dict = assign_anomaly(out)
     pprint.pprint(set(anomaly_dict.values()))


    #=====================================================





