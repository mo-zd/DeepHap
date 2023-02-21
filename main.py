############################ main.py ##################################
#   This is a python package for analysing hap data
#   Author: Mohsen Yazdani
#   Email: yazdani.m@ut.ac.ir
#######################################################################

import copy
import sys

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib
matplotlib.use('Qt5Agg')
from PyQt5.QtWidgets import QFileDialog
from PyQt5 import QtCore, QtWidgets

from plot import *
from util import *
from gui.haplo import *


class MainDialog(QtWidgets.QDialog):

    def __init__(self, *args, **kwargs):
        self.files = {}
        self.figs = {}
        self.window_size = 700
        super(MainDialog, self).__init__(*args, **kwargs)
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        self.ui.m1_button.clicked.connect(self.get_m1_file)
        self.ui.m2_button.clicked.connect(self.get_m2_file)
        self.ui.p1_button.clicked.connect(self.get_p1_file)
        self.ui.p2_button.clicked.connect(self.get_p2_file)
        self.ui.logr_button.clicked.connect(self.get_logr_file)
        self.ui.analyze_button.clicked.connect(self.analyze)
        self.show()

    def get_m1_file(self):
        fname, _ = QFileDialog.getOpenFileName(self, 'Open File')
        self.files['m1'] = fname

    def get_m2_file(self):
        fname, _ = QFileDialog.getOpenFileName(self, 'Open File')
        self.files['m2'] = fname

    def get_p1_file(self):
        fname, _ = QFileDialog.getOpenFileName(self, 'Open File')
        self.files['p1'] = fname

    def get_p2_file(self):
        fname, _ = QFileDialog.getOpenFileName(self, 'Open File')
        self.files['p2'] = fname

    def get_logr_file(self):
        fname, _ = QFileDialog.getOpenFileName(self, 'Open File')
        self.files['logr'] = fname

    def read_window_size(self):
        self.window_size = int(self.ui.window_size_edit.text())

    def show_plot(self):
        self.figs['m1'].setParent(self.ui.frame_m1)
        self.figs['m1'].show()
        self.figs['m2'].setParent(self.ui.frame_m2)
        self.figs['m2'].show()
        self.figs['p1'].setParent(self.ui.frame_p1)
        self.figs['p1'].show()
        self.figs['p2'].setParent(self.ui.frame_p2)
        self.figs['p2'].show()
        self.figs['logr'].setParent(self.ui.frame_logr)
        self.figs['logr'].show()

    def analyze(self):
        self.read_window_size()
        window_size = self.window_size
        m1_data = read_baf_file(self.files['m1'])
        m1_info = process_baf_data(m1_data, window_size=window_size)
        m2_data = read_baf_file(self.files['m2'])
        m2_info = process_baf_data(m2_data, window_size=window_size)
        p1_data = read_baf_file(self.files['p1'])
        p1_info = process_baf_data(p1_data, window_size=window_size)
        p2_data = read_baf_file(self.files['p2'])
        p2_info = process_baf_data(p2_data, window_size=window_size)
        logr_data = read_logr_file(self.files['logr'])
        logr_info = process_logr_data(logr_data, window_size=window_size)

        out = generate_hap_signatures(
            m1_info, m2_info, p1_info, p2_info, logr_info)
        anomaly_dict = assign_anomaly(out)
        pprint.pprint(set(anomaly_dict.values()))
        to_print = '\n'.join(list(set(anomaly_dict.values())))
        self.ui.output_label.setText(to_print)
        self.show()
        # print('timing:', time.time() - start)

        self.figs['m1'] = MplCanvas()
        self.figs['m2'] = MplCanvas()
        self.figs['p1'] = MplCanvas()
        self.figs['p2'] = MplCanvas()
        self.figs['logr'] = MplCanvas()

        plot_logr_gui(logr_data, logr_info, self.figs['logr'])
        plot_baf_gui(m1_data, m1_info, self.figs['m1'])
        plot_baf_gui(m2_data, m2_info, self.figs['m2'])
        plot_baf_gui(p1_data, p1_info, self.figs['p1'])
        plot_baf_gui(p2_data, p2_info, self.figs['p2'])

        # self.ui.plot_type_comboBox.currentText().lower()
        # self.figs['m1'].setParent(self.ui.frame)
        # self.figs['m1'].show()
        self.show_plot()


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    w = MainDialog()
    app.exec_()
