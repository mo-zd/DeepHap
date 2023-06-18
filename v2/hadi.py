from collections import Counter
import time

import numpy as np
import matplotlib.pyplot as plt

from settings import BAF_ASSIGNMENT_DICT
from settings import LOGR_ASSIGNMENT_DICT
from settings import CENTROMER_DICT
from settings import ANOMALY_DICT
from read_write import *
from exceptions import *


class Denoiser(Base):
    def __init__(self, file_reader, window_size=600):
        super().__init__()
        self.fr = file_reader
        self.input_data = file_reader.data
        self.denoised = dict(M1=None, M2=None, P1=None, P2=None, Logr=None)
        self.window_size = window_size

    def _denoise(self, values, allowed_values):
        starts = [0 for _ in range(int(self.window_size / 2))] + list(range(0, len(values) - (int(self.window_size / 2))))
        ends = list(range(int(self.window_size / 2), len(values))) + [(len(values) - 1) for _ in range(int(self.window_size / 2))]
        windows = list(zip(starts, ends))[:len(values)]
        out = []
        for window in windows:
            counts = Counter(values[window[0]: window[1]])
            val = allowed_values[np.argmax([counts[i] for i in allowed_values])]
            out.append(val)
        return out

    def _denoise_single_baf(self, df):
        col_name = df.columns[-1]
        global BAF_ASSIGNMENT_DICT
        df_list = []
        chrs = self.get_chrs(df)
        df_list = self.split(df)
        df_list_denoised = []
        for chr in chrs:
            df_chr = df_list[self.order_dict[chr]].copy()
            baf_values = df_chr[col_name].tolist()
            if baf_values[0] == -1:
                df_chr['denoised'] = baf_values
                df_list_denoised.append(df_chr)
            else:
                baf_values_denoised = [BAF_ASSIGNMENT_DICT[int(baf * 12.)] for baf in baf_values]
                allowed_values = sorted(list(set(baf_values_denoised)))
                if self.window_size > 0:
                    baf_values_denoised = self._denoise(baf_values_denoised, allowed_values)
                df_chr['denoised'] = baf_values_denoised
                df_list_denoised.append(df_chr)
        df = self.concat(df_list_denoised)
        return df
        
    def denoise_baf(self):
        for label in self.input_data:
            if label in ('M1', 'M2', 'P1', 'P2'):
                df = self.input_data[label]
                df = self._denoise_single_baf(df)
                self.denoised[label] = df
            
    def denoise_logr(self):
        global LOGR_ASSIGNMENT_DICT
        df = self.input_data['Logr']
        col_name = df.columns[-1]
        df_list = self.split(df)
        chrs = self.get_chrs(df)
        df_list_denoised = []
        for chr in chrs:
            df_chr = df_list[self.order_dict[chr]].copy()
            logr_values = df_chr[col_name].tolist()
            logr_values = np.clip(np.array(logr_values), -2, 0).tolist()
            if logr_values[0] == -1:
                df_chr['denoised'] = logr_values
                df_list_denoised.append(df_chr)
            else:
                logr_values_denoised = [LOGR_ASSIGNMENT_DICT[int((logr + 2.) * 2.)] for logr in logr_values]
                allowed_values = sorted(list(set(logr_values_denoised)))
                if self.window_size > 0:
                    logr_values_denoised = self._denoise(logr_values_denoised, allowed_values)
                df_chr['denoised'] = logr_values_denoised
                df_list_denoised.append(df_chr)
        df = self.concat(df_list_denoised)
        self.denoised['Logr'] = df
      
class Assigner(Base):
    def __init__(self, denoiser):
        super().__init__()
        self.denoiser = denoiser
        self.denoised = denoiser.denoised

    def _nearest(self, pivot, query):
        out = []
        lb, ub = -1000, query[0]
        mid = (lb + ub) / 2.
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
                out.append((el, lb))
            else:
                out.append((el, ub))
        return out

    def _align_chr(self, chr):
        chr_index = self.denoiser.order_dict[chr]
        denoised_chr = dict([(el, self.denoiser.split(self.denoised[el])[chr_index]) for el in ('M1', 'M2', 'P1', 'P2', 'Logr')])
        pos_chr = dict([(el, denoised_chr[el]['Position'].tolist()) for el in ('M1', 'M2', 'P1', 'P2', 'Logr')])
        all_aligns = []
        for pivot in ('M1', 'M2', 'P1', 'P2', 'Logr'):
            pivot_aligns = []
            for query in ('M1', 'M2', 'P1', 'P2', 'Logr'):
                _, q = zip(*self._nearest(pos_chr[pivot], pos_chr[query]))
                pivot_aligns.append(q)
            pivot_aligns = list(zip(*pivot_aligns))
        all_aligns.extend(pivot_aligns)

        #removing redundant alignmnents and sorting
        all_aligns = list(set(all_aligns))
        all_aligns = list(filter(lambda el: min(el) > 0, all_aligns))
        all_aligns = sorted(all_aligns)
        return all_aligns
    
    def align(self):
        alignment = dict()
        chrs_dict = dict()
        chrs = set(self.get_chrs(self.denoised['M1']))
        for el in ('M1', 'M2', 'P1', 'P2', 'Logr'):
            chrs_dict[el] = self.get_chrs(self.denoised[el])
            chrs = chrs.intersection(chrs_dict[el])
        chrs = list(chrs)
        chrs = sorted(chrs, key=lambda el: self.order_dict[el])
        if 'Y' in chrs:
            chrs.remove('Y')
        for chr in chrs:
            alignment[chr] = self._align_chr(chr)
        self.alignment = alignment

    def assign_anomaly(self):
        anomalies = dict()
        signatures = dict()
        for chr in self.alignment.keys():
            signatures_chr = []
            denoised_chr = dict()
            pos_dict_chr = dict()
            for el in ('M1', 'M2', 'P1', 'P2', 'Logr'):
                denoised_chr[el] = self.split(self.denoised[el])[self.order_dict[chr]]
                pos = list(denoised_chr[el]['Position'].tolist())
                val = list(denoised_chr[el]['denoised'].tolist())
                pos_dict_chr[el] = dict(zip(pos, val))
            alignment_chr = self.alignment[chr]
            anomalies_chr = []
            for el in alignment_chr:
                try:
                    signature = (pos_dict_chr['M1'][el[0]],
                                 pos_dict_chr['M2'][el[1]],
                                 pos_dict_chr['P1'][el[2]],
                                 pos_dict_chr['P2'][el[3]],
                                 pos_dict_chr['Logr'][el[4]])
                    signatures_chr.append(signature)
                    anomalies_chr.append((el, ANOMALY_DICT[signature]))
                except:
                    pass
            anomalies[chr] = anomalies_chr
            signatures[chr] = signatures_chr
        self.anomalies = anomalies
        self.signatures = signatures
        return anomalies
              
