import pandas as pd
import os.path as osp


def read_baf_file(baf_file):
    df = pd.read_csv(baf_file, header=None, sep='\t')
    baf_data = {}
    for ind, chr_name, loc, baf in df.values:
        baf_data[ind] = (chr_name, loc, baf)
    return baf_data


def read_logr_file(logr_file):
    df = pd.read_csv(logr_file, header=None, sep='\t')
    logr_data = {}
    for ind, chr_name, loc, cval in df.values:
        logr_data[ind] = (chr_name, loc, cval)
    return logr_data


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
    f_name = osp.join(osp.dirname(__file__), 'files','Tbl_M1_Seg.txt')
    extract_sample_files(f_name)


