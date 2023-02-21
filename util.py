############################ util.py ##################################
#   This is a python package for analysing hap data
#   Author: Mohsen Yazdani
#   Email: yazdani.m@ut.ac.ir
#######################################################################

import copy
import time
import os.path as osp
from collections import Counter
import pprint
import random

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from setting import BAF_ASSIGNMENT_DICT
from setting import LOGR_ASSIGNMENT_DICT
from setting import CENTROMER_DICT
from setting import ANOMALY_DICT
from plot import plot_baf
from plot import plot_logr
from read_write import read_baf_file_all
from read_write import read_logr_file_all


#TODO: work on denoise function make it better
def denoise(array, allowed_values, window_size=600):
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
    value = list(list(baf_data.values())[0].values())[0][2]
    if value == -1:
        return None
    global BAF_ASSIGNMENT_DICT
    baf_info = dict()
    for chr_name, _baf_data in baf_data.items():
        _baf_info = dict()
        for ind, (chr_name, loc, baf) in _baf_data.items():
            _baf_info[ind] = (chr_name, loc, BAF_ASSIGNMENT_DICT[int(baf * 12.)])
        if window_size:
            out = sorted(list(_baf_info.items()), key=lambda x: x[1][1])
            keys, vals = list(zip(*out))
            chr_names, locs, bafs = zip(*vals)
            allowed_values = sorted(tuple(set(bafs)))
            bafs = denoise(bafs, allowed_values, window_size)
            _baf_info = dict(zip(keys, list(zip(chr_names, locs, bafs))))
        baf_info[chr_name] = _baf_info
    return baf_info


def process_logr_data(logr_data, window_size=0):
    global LOGR_ASSIGNMENT_DICT
    logr_info = dict()
    for chr_name, _logr_data in logr_data.items():
        _logr_info = dict()
        for ind, (chr_name, loc, logr) in _logr_data.items():
            _logr_info[ind] = (chr_name, loc, LOGR_ASSIGNMENT_DICT[int((logr + 2.) * 2.)])
        if window_size:
            out = sorted(list(_logr_info.items()), key=lambda x: x[1][1])
            keys, vals = list(zip(*out))
            chr_names, locs, logrs = zip(*vals)
            allowed_values = sorted(tuple(set(logrs)))
            logrs = denoise(logrs, allowed_values, window_size)
            _logr_info = dict(zip(keys, list(zip(chr_names, locs, logrs))))
        logr_info[chr_name] = _logr_info
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

#TODO: change below function to be compatible with new *_info values
# def generate_neighbor_locs(m1_info, m2_info, p1_info, p2_info, logr_info):
#     m1 = [(k, v[1]) for k, v in m1_info.items()]
#     m2 = [(k, v[1]) for k, v in m2_info.items()]
#     p1 = [(k, v[1]) for k, v in p1_info.items()]
#     p2 = [(k, v[1]) for k, v in p2_info.items()]
#     logr = [(k, v[1]) for k, v in logr_info.items()]

#     m1_list = [el[0] for el in m1]
#     m2_list = [el[0] for el in m2]
#     p1_list = [el[0] for el in p1]
#     p2_list = [el[0] for el in p2]
#     logr_list = [el[0] for el in logr]

#     set1 = set(zip(m1_list, nearest(m2, m1), nearest(p1, m1), nearest(p2, m1), nearest(logr, m1)))
#     set2 = set(zip(nearest(m1, m2), m2_list, nearest(p1, m2), nearest(p2, m2), nearest(logr, m2)))
#     set3 = set(zip(nearest(m1, p1), nearest(m2, p1), p1_list, nearest(p2, p1), nearest(logr, p1)))
#     set4 = set(zip(nearest(m1, p2), nearest(m2, p2), nearest(p1, p2), p2_list, nearest(logr, p2)))
#     set5 = set(zip(nearest(m1, logr), nearest(m2, logr), nearest(p1, logr), nearest(p2, logr), logr_list))
#     out = set1.union(set2, set3, set4, set5)
#     return out

def generate_neighbor_locs(*info_list):
    _info_list = [el for el in info_list if el]
    neighbor_locs = dict()
    chrs = [str(el) for el in range(1, 23)]
    chrs.append('X')
    for chr in chrs:
        d1s = []
        d1_lists = []
        for info in _info_list:
            d1 = [(k, v[1]) for k, v in info[chr].items()]
            d1_list = [el[0] for el in d1]
            d1s.append(d1)
            d1_lists.append(d1_list)
        sets = []
        for i in range(len(d1s)):
            temp = []
            temp.append(d1_lists[i])
            indices = list(range(len(d1s)))
            indices.remove(i)
            for j in indices:
                temp.append(nearest(d1s[j], d1s[i]))
            sets.append(set(zip(*temp)))
        out = sets[0].union(*sets[1:])
        neighbor_locs[chr] = out   
    return neighbor_locs  



#TODO: change below function to be compatible with new *_info values
# def generate_hap_signatures(m1_info, m2_info, p1_info, p2_info, logr_info):
#     out = generate_neighbor_locs(m1_info, m2_info, p1_info, p2_info, logr_info)
#     out = {(e1, e2, e3, e4, e5): (m1_info[e1][2], m2_info[e2][2], p1_info[e3][2], p2_info[e4][2], logr_info[e5][2]) for
#      e1, e2, e3, e4, e5 in out}
#     return out

def generate_hap_signatures(info_list):
    out = generate_neighbor_locs(*info_list)
    chrs = [str(el) for el in range(1, 23)]
    chrs.append('X')
    result = dict()
    for chr in chrs:
        out[chr] = {(e1, e2, e3, e4, e5): [info_list] for e1, e2, e3, e4, e5 in out[chr]}
    return out

#TODO: change below function to be compatible with new hap_signatures
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
    start = time.time()
    window_size = 600
    logr_file = osp.join(osp.dirname(__file__), 'BAF',
            'Mother.Embryo_3794p1_eA2_2.5.LogR.csv')
    m1_file = osp.join(osp.dirname(__file__), 'BAF',
            'Mother.Embryo_3794p1_eA2_2.5.M1.csv')
    m2_file = osp.join(osp.dirname(__file__), 'BAF',
            'Mother.Embryo_3794p1_eA2_2.5.M2.csv')
    p1_file = osp.join(osp.dirname(__file__), 'BAF',
            'Mother.Embryo_3794p1_eA2_2.5.P1.csv')
    p2_file = osp.join(osp.dirname(__file__), 'BAF',
            'Mother.Embryo_3794p1_eA2_2.5.P2.csv')
    files = {'m1': m1_file,
             'm2': m2_file,
             'p1': p1_file,
             'p2': p2_file,
             'logr': logr_file}
    
    data = dict()
    for k in ('m1', 'm2', 'p1', 'p2'):
        data[k] = read_baf_file_all(files[k])
    data["logr"] = read_logr_file_all(files['logr'])

    info = dict()
    for k in ("m1", "m2", "p1", "p2"):
        info[k] = process_baf_data(data[k], window_size=window_size)
    info["logr"] = process_logr_data(data[k], window_size=window_size)

    out = generate_neighbor_locs(*list(info.values()))
    print(len(out))
    print(len(out['1']))

    #find missing parent

    # missing = set()
    # for k in ('m1', 'm2', 'p1', 'p2'):
    #     if not info[k]:
    #         missing.add(k[0])

    # if missing:
    #     missing_parent = missing.pop()
    #     avai_parent = ({'m', 'p'} - {missing_parent}).pop()
    #     for k in ('1', '2'):
    #         info[missing_parent + k] = data[missing_parent + k]
    #         avai_info = info[avai_parent + k] 
    #         for chr_name, baf_info in avai_info.items():
    #             temp = list(filter(lambda el: CENTROMER_DICT[chr_name][0]<el[1][1]<CENTROMER_DICT[chr_name][1], list(baf_info.items())))
    #             if len(temp):
    #                 temp = [el[1][2] for el in temp]
    #                 value = max(set(temp), key=temp.count)
    #                 print(chr_name, value, len(temp))
    #                 for k1, v in info[missing_parent + k][chr_name].items():
    #                     info[missing_parent + k][chr_name][k1] = [v[0], v[1], value]
    #             else:
    #                 print(chr_name, "not found")   


    # print("m1_info:", m1_info)
    # print("m2_info:", m2_info)
    # print("p1_info:", p1_info)
    # print("p2_info:", p2_info)
    # print("logr_info:", logr_info)


    #  plot_logr(logr_data, logr_info)
    #  plot_baf(m1_data, m1_info)
    #  plot_baf(m2_data, m2_info)
    #  plot_baf(p1_data, p1_info)
    #  plot_baf(p2_data, p2_info)


    #  out = generate_hap_signatures(m1_info, m2_info, p1_info, p2_info, logr_info)
    #  print('timing:', time.time() - start)
    #  print('number of elements of out is', len(out))
    #  anomaly_dict = assign_anomaly(out)
    #  pprint.pprint(set(anomaly_dict.values()))
