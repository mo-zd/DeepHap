import os.path as osp
import pandas as pd
import pickle

from exceptions import *


class Base(object):
    def __init__(self):
        self.order_dict = dict(zip([str(i) for i in range(1, 23)] + ['X', 'Y'], range(0, 24)))

    def get_chrs(self, df):
        
        # Add NaN handling
        df = df.dropna(subset=['Chr'])
        
        chrs = list(set(df['Chr'].tolist()))
        chrs = sorted(chrs, key=lambda el: self.order_dict[el])
        self.chrs = chrs
        return chrs
   

    def concat(self, df_list):
        df = pd.concat(df_list, ignore_index=True)
        return df
    
    def split(self, df):
        chrs = self.get_chrs(df)
        df_list = []
        for chr in chrs:
            df_chr = df[df['Chr'] == chr]
            df_list.append(df_chr)
        return df_list

    def save(self, file_path):
        with open(file_path, "wb") as f:
            pickle.dump(self, f)


class Reader(Base):
    def __init__(self):
        self.data = dict(M1=None, M2=None, P1=None, P2=None, Logr=None)

    def read_baf(self, baf_file, label):
        if not(label in ('M1', 'M2', 'P1', 'P2')):
            raise FileLabelError
        dtype = dict()
        dtype["Chr"] = str
        df = pd.read_csv(baf_file, dtype=dtype)
        columns = list(df.keys())
        col_name = list(filter(lambda el: el[0] == 'E', columns))[0]
        df = df[["Chr", "Position", col_name]]
        self.data[label] = df

    def read_logr(self, logr_file):
        dtype = {"Chr": str}
        df = pd.read_csv(logr_file, dtype=dtype)
        columns = list(df.keys())
        col_name = list(filter(lambda el: el[0] == 'E', columns))[0]
        df = df[["Chr", "Position", col_name]]
        self.data['Logr'] = df

    def extract_sample_files(self, f_name):
        df = pd.read_csv(f_name, sep='\t')
        for sample in df.keys():
            if sample[0] == 'E':
                df1 = df[['Name', 'Chr', 'Position', sample]]
                df1.to_csv(
                    osp.join(
                        osp.dirname(f_name),
                        osp.basename(f_name)[:-4] + f'_{sample}.csv'
                    ), sep=',', index=None)


def load(file_path):
    with open(file_path, "rb") as f:
        out = pickle.load(f)
    return out
