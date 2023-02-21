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