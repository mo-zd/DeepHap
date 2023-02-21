######################### read_write.py ###############################
#   This is a python package for analysing hap data
#   Author: Mohsen Yazdani
#   Email: yazdani.m@ut.ac.ir
#######################################################################

import json

import pandas as pd
import os.path as osp
from itertools import groupby


# def read_baf_file(baf_file):
#     df = pd.read_csv(baf_file, header=None, sep='\t')
#     baf_data = {}
#     for ind, chr_name, loc, baf in df.values:
#         baf_data[ind] = (chr_name, loc, baf)
#     return baf_data


# def read_logr_file(logr_file):
#     df = pd.read_csv(logr_file, header=None, sep='\t')
#     logr_data = {}
#     for ind, chr_name, loc, cval in df.values:
#         logr_data[ind] = (chr_name, loc, cval)
#     return logr_data


def read_baf_file_all(baf_file):
    '''This is function for reading baf file'''
    dtype = dict()
    dtype["Chr"] = str
    df = pd.read_csv(baf_file, dtype=dtype)
    baf_data = {}
    columns = list(df.keys())
    col_name = list(filter(lambda el: el[0] == 'E', columns))[0]
    df1 = df[["Name", "Chr", "Position", col_name]]
    for ind, chr_name, loc, baf in df1.values:
        baf_data[ind] = (str(chr_name), int(loc), baf)
    baf_data = list(baf_data.items())
    baf_data = sorted(baf_data, key=lambda el: el[1][0])
    groups = groupby(baf_data, key=lambda el: el[1][0])
    out = dict()
    for k, g in groups:
        out[k] = dict(list(g))
    return out


def read_logr_file_all(logr_file):
    dtype = dict()
    dtype["Chr"] = str
    df = pd.read_csv(logr_file, dtype=dtype)
    logr_data = {}
    columns = list(df.keys())
    col_name = list(filter(lambda el: el[0] == 'E', columns))[0]
    df1 = df[["Name", "Chr", "Position", col_name]]
    logr_file = {}
    for ind, chr_name, loc, baf in df1.values:
        logr_data[ind] = (str(chr_name), int(loc), baf)
    logr_data = list(logr_data.items())
    logr_data = sorted(logr_data, key=lambda el: el[1][0])
    groups = groupby(logr_data, key=lambda el: el[1][0])
    out = dict()
    for k, g in groups:
        out[k] = dict(list(g))
    return out


def extract_sample_files(f_name):
    df = pd.read_csv(f_name, sep='\t')
    for sample in df.keys():
        if sample[0] == 'E':
            df1 = df[['Name', 'Chr', 'Position', sample]]
            df1.to_csv(
                osp.join(
                    osp.dirname(f_name),
                    osp.basename(f_name)[:-4] + f'_{sample}.txt'
                ), sep='\t', index=None)


if __name__ == '__main__':
    f_name = osp.join(osp.dirname(__file__), 'BAF',
                      'Mother.Embryo_3794p1_eA2_2.5.M1.csv')
    baf_data = read_baf_file_all(f_name)
    print(baf_data['X'])

