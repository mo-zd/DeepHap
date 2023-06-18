import argparse
import os
import sys
import glob
import pickle

from hadi import *

VERSION = '0.0.1'
parser = argparse.ArgumentParser(description='Finding Anomalies')
parser.version = VERSION
parser.add_argument('i',
                       metavar='input_folder',
                       type=str,
                       help='the path to input folder containing M1, M2, P1, P2, LogR csv files')

parser.add_argument('o', metavar='output_folder', type=str, help="the path to output folder")
parser.add_argument('-w', '--window_size', metavar='window_size', type=int, help="window size for denoising operator", default=600)
parser.add_argument('-v', '--version', action='version', help="package version")

args = parser.parse_args()

input_path = args.i
if not os.path.isdir(input_path):
    print('The path specified does not exist')
    sys.exit()

output_path = args.o
if not os.path.isdir(output_path):
    os.mkdir(output_path)

# running
files = glob.glob(os.path.join(input_path, '*.csv'))

files_dict = dict()
for file_type in ('m1', 'm2', 'p1', 'p2', 'logr'):
    candid_files = list(filter(lambda f: f.lower().endswith(f"{file_type}.csv"), files))
    if len(candid_files) == 1:
        files_dict[file_type.upper()] = candid_files[0]
    else:
        print(f'ERROR: {file_type.upper()} file does not exist. Please add this file to input folder.')
        sys.exit()

start = time.time()

fr = Reader()
for file_type in ('M1', 'M2', 'P1', 'P2'):
    fr.read_baf(files_dict[file_type], file_type)
fr.read_logr(files_dict['LOGR'])

den = Denoiser(fr, window_size=args.window_size)
den.denoise_baf()
den.denoise_logr()

assigner = Assigner(den)
alignment = assigner.align()
anomalies = assigner.assign_anomaly()

with open(os.path.join(output_path, "model"), "wb") as f:
    pickle.dump({'denoiser': den, 'assigner': assigner}, f)

# saving results
with open(os.path.join(output_path, "anomalies.txt"), 'w') as f:
    for k, v in assigner.anomalies.items():
        try:
            _, anom = zip(*v)
            anom = sorted(list(set(anom)))
            f.write(f"{k}:\t\t{', '.join(anom)}\n")
        except:
            pass

print(f"finished in {round(time.time() - start)} seconds!")

