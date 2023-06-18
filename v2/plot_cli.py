import argparse
import os
import sys
import pickle

from plots import * 

VERSION = '0.0.1'
parser = argparse.ArgumentParser(description="plot chromosomes")
parser.version = VERSION
parser.add_argument('o', metavar='output_folder', type=str, help="the path to output folder")
parser.add_argument('c', metavar='chromosome', type=str, help="chromosome name")

args = parser.parse_args()
output_path = args.o
if not (os.path.isdir(output_path)):
    print('Errot: The path specified does not exist')
    sys.exit()

if not(os.path.isfile(os.path.join(output_path, "model"))):
    print('ERROR: The model file does not exist')
    sys.exit()

with open(os.path.join(output_path, 'model'), "rb") as f:
    model = pickle.load(f)
den = model['denoiser']
assigner = model['assigner']
chr = args.c.upper()
plot(den, chr)
plt.show()

    
