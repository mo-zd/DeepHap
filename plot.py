############################ plot.py ##################################
#   This is a python package for analysing hap data
#   Author: Mohsen Yazdani
#   Email: yazdani.m@ut.ac.ir
#######################################################################

import json
import os.path as osp

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None):
        fig = Figure()
        self.ax1 = fig.add_subplot(211)
        self.ax2 = fig.add_subplot(212)
        super(MplCanvas, self).__init__(fig)

def plot_baf_gui(baf_data, baf_info, mplCanvas):
    ax1 = mplCanvas.ax1
    baf_info_vals = [val[2] for val in baf_info.values()]
    ax1.scatter(range(1, len(baf_info_vals) + 1), baf_info_vals, c='b', s=0.2)
    plt.title('BAF Information')
    plt.yticks([0, 1./3, 1./2, 2./3, 1.])
    plt.grid(True, alpha=0.3)

    ax2 = mplCanvas.ax2
    baf_data_vals = [val[2] for val in baf_data.values()]
    ax2.scatter(range(1, len(baf_data_vals) + 1), baf_data_vals, c='r', s=0.2)
    plt.title('BAF Data')
    plt.yticks([0, 1./3, 1./2, 2./3, 1.])
    plt.grid(True, alpha=0.3)

def plot_logr_gui(logr_data, logr_info, mplCanvas):
    ax1 = mplCanvas.ax1
    logr_info_vals = [val[2] for val in logr_info.values()]
    ax1.scatter(range(1, len(logr_info_vals) + 1), logr_info_vals, c='b', s=0.2)
    plt.title('Copy Number Information')
    plt.yticks([0, 1, 2, 3, 4])
    plt.grid(True, alpha=0.3)

    ax2 = mplCanvas.ax2
    logr_data_vals = [val[2] for val in logr_data.values()]
    ax2.scatter(range(1, len(logr_data_vals) + 1), logr_data_vals, c='r', s=0.2)
    plt.title('LogR Data')
    plt.yticks([-3, -2, -1, 0, 1])
    plt.grid(True, alpha=0.3)

def plot_baf(baf_data, baf_info):
    fig = plt.figure()

    ax1 = plt.subplot(211)
    baf_info_vals = [val[2] for val in baf_info.values()]
    ax1.scatter(range(1, len(baf_info_vals) + 1), baf_info_vals, c='b', s=0.2)
    plt.title('BAF Information')
    plt.yticks([0, 1./3, 1./2, 2./3, 1.])
    plt.grid(True, alpha=0.3)

    ax2 = plt.subplot(212)
    baf_data_vals = [val[2] for val in baf_data.values()]
    ax2.scatter(range(1, len(baf_data_vals) + 1), baf_data_vals, c='r', s=0.2)
    plt.title('BAF Data')
    plt.yticks([0, 1./3, 1./2, 2./3, 1.])
    plt.grid(True, alpha=0.3)
    return fig

def plot_logr(logr_data, logr_info):
    fig = plt.figure()

    ax1 = plt.subplot(211)
    logr_info_vals = [val[2] for val in logr_info.values()]
    ax1.scatter(range(1, len(logr_info_vals) + 1), logr_info_vals, c='b', s=0.2)
    plt.title('Copy Number Information')
    plt.yticks([0, 1, 2, 3, 4])
    plt.grid(True, alpha=0.3)

    ax2 = plt.subplot(212)
    logr_data_vals = [val[2] for val in logr_data.values()]
    ax2.scatter(range(1, len(logr_data_vals) + 1), logr_data_vals, c='r', s=0.2)
    plt.title('LogR Data')
    plt.yticks([-3, -2, -1, 0, 1])
    plt.grid(True, alpha=0.3)
    return fig

#def plot_genome_logr()

def plot_baf_1(baf_data1, baf_data2, root_dir):
    data1 = {}
    chrs = [str(el) for el in range(1, 23)] + ['X']
    for chr in chrs:
        data1[chr] = []
        for k, v in baf_data1[chr].items():
            el = {
                "x": v[1],
                "y": v[2],
            }
            data1[chr].append(el)
    
    data2 = {}
    for chr in chrs:
        data2[chr] = []
        for k, v in baf_data2[chr].items():
            el = {
                "x": v[1],
                "y": v[2],
            }
            data2[chr].append(el)
    data = [data1, data2]
    path = osp.join(root_dir, "js", "data.js")
    data = 'var data = ' + f" \'{json.dumps(data)}\'"
    with open(path, "w") as f:
        f.write(data)


if __name__ == '__main__':
    root_dir = '/media/mohsen/F4CE1B6ECE1B27FE/projects/yazdi_thesis/haplotyping'
    from read_write import read_baf_file_all
    baf_data1 = read_baf_file_all(osp.join(root_dir, "BAF", "Mother.Embryo_3794p1_eA2_2.5.M1.csv"))
    baf_data2 = read_baf_file_all(osp.join(root_dir, "BAF", "Mother.Embryo_3794p1_eA2_2.5.M2.csv"))
    plot_baf_1(baf_data1, baf_data2,  root_dir)

