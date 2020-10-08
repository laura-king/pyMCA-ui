'''
File responsonsible for reading and processing data for SSRL_MCA pydm screen.

Original:
Refactored: Laura King
09/2020
'''
import math
import argparse
import logging
import numpy as np
import scipy.signal as signal

from qtpy.QtCore import Qt
from qtpy import QtWidgets
from pydm import Display
from pydm.widgets.channel import PyDMChannel
from scipy.optimize import curve_fit
from operator import itemgetter
from os import path


def build_dic(cli_args):
    emission = {}
    element = {}
    filename = "XRay_Emission_Lines.txt"

    if cli_args.filename:
        filename = cli_args.filename

    with open(filename) as file:
        for line in file:
            if (len(line) <= 1):
                continue
            line_data = line.split()
            try:
                int(line_data[0])
            except BaseException:
                continue

            # First two always symbols
            symbol = line_data[1] + line_data[0]

            # Convert numbers to floats, set all dashes to and pad with 0's
            energy = []
            for energy_val in line_data[2:]:
                try:
                    energy.append(float(energy_val))
                except BaseException:
                    energy.append(0.)
            energy = energy + [0.] * (9 - len(energy))

            # Set vars to data
            element_d = {}
            elem_dict = {
                'Ka': [energy[0], energy[1]],
                'Kb': [energy[2]],
                'La': [energy[3], energy[4]],
                'Lb': [energy[5], energy[6]],
                'Lg': [energy[7]],
                'Ma': [energy[8]]
            }
            # Filter out 0's, assign data to dictionaries
            for elem, value_list in elem_dict.items():
                full_symbol = symbol + "-" + elem
                if len(value_list) == 2:
                    if value_list[0] != 0 and value_list[1] != 0:
                        element_d, emission = combine(
                            (full_symbol), value_list, element_d, emission)
                        continue
                    else:
                        element_d, emission = add_to_dicts(
                            (full_symbol + "2"), value_list[1], element_d, emission)
                element_d, emission = add_to_dicts(
                    (full_symbol + "1"), value_list[0], element_d, emission)
            element[symbol] = element_d

        energy_i = []
        for key in emission.keys():
            keys = key.split("-")
            energy_i.append([emission[key], keys[0], keys[1]])

        sort_energy = sorted(energy_i, key=itemgetter(0))
        return sort_energy, element


def add_to_dicts(symbol, value, element_d, emission):
    """
    Adds the value-symbol pair to the element and emission dictionaries
    """
    if value != 0:
        name = symbol.split('-')[-1]
        element_d[name] = value
        emission[symbol] = value
    return element_d, emission


def combine(symbol, value_list, element_d, emission):
    """
    Adds the correct symbol-value pairs for the double variables to the
    element and emission dictionaries
    """
    if ((value_list[0] - value_list[1]) > -
            30) and ((value_list[0] - value_list[1]) < 30):
        element_d, emission = add_to_dicts(
            symbol, ((value_list[0] + value_list[1]) / 2.), element_d, emission)
    else:
        element_d, emission = add_to_dicts(
            symbol + "1", value_list[0], element_d, emission)
        element_d, emission = add_to_dicts(
            symbol + "2", value_list[1], element_d, emission)
    return element_d, emission


def gaussian(x, amplitude, mean, sigma):
    return amplitude * np.exp(-((x - mean) / sigma)**2 / 2.)


def gaussian2(x, ampl1, mean1, sigma1, ampl2, mean2, sigma2):
    return ampl1 * np.exp(-((x - mean1) / sigma1)**2 / 2.) +                         \
        ampl2 * np.exp(-((x - mean2) / sigma2)**2 / 2.)


def cauchy(x, amplitude, mean, gamma):
    return amplitude * gamma / ((x - mean) * (x - mean) + gamma * gamma)


class MCADisplay(Display):
    def __init__(self, parent=None, args=None, macros=None):
        super(
            MCADisplay,
            self).__init__(
            parent=parent,
            args=args,
            macros=macros)
        # Debug Logger
        self.logger = logging.getLogger('mca_logger')
        self.separator = "\n" + ("-" * 20) + "\n"

        self.num_ROI = 9
        self.ROI = []
        self.start = []
        self.end = []
        self.counts = []
        self.lines = []
        self.set_ROI_widgets()
        cli_args = self.parse_args(args)

        self.energy, self.element = build_dic(cli_args)

        self.waveform.plotItem.scene().sigMouseMoved.connect(self.mouse_moved)
        self.waveform.setXLabels(["Energy (eV)"])
        self.waveform.setYLabels(["Count"])

        # Add Channels
        self.waveform.addChannel(None, None, name="Full", color="white")
        color_list = ["red", "green", "blue"]
        for wave in range(self.num_ROI):
            name = f"ROI{wave+1}"
            color = color_list[wave % len(color_list)]
            self.waveform.addChannel(
                None, None, name=name, color=color, lineWidth=2)

        # TODO: Is 18 just double the number of ROI's?
        for wave in range(18):
            name = f"Line{wave+1:02d}"
            self.waveform.addChannel(
                None,
                None,
                name=name,
                color="white",
                lineWidth=2,
                lineStyle=Qt.DashLine)

        self.curve = self.waveform._curves[0]
        self.croi = self.waveform._curves[1:10]
        self.line = self.waveform._curves[10:28]

        if (macros is not None) and ("FIT" in macros):
            if (macros["FIT"].lower() == "cauchy"):
                self.fitc = "Cauchy"
            else:
                self.fitc = "Gaussian"
        else:
            self.fitc = "Gaussian"

        if (macros is not None) and ("DEVICE" in macros):
            self.dataSource.addItem("Live EPICS")

            # TODO: Other file uses macros["DEVICE"]+":RAW:ArrayData"
            self.epics = PyDMChannel(address="ca://" +
                                macros["DEVICE"] + ":ARR1:ArrayData",
                                value_slot=self.live_data)
            self.epics.connect()

            self.show_exposure()
        else:
            self.show_mca()

        self.dataSource.addItem("Playback")
        self.dataSource.setCurrentIndex(0)
        self.dataSource.currentIndexChanged.connect(self.change_source)

        self.connectButton.clicked.connect(self.connect_data)
        self.openFile.clicked.connect(self.open_file)
        self.previousMCA.clicked.connect(self.previous_mca)
        self.nextMCA.clicked.connect(self.next_mca)
        self.fullView.clicked.connect(self.full_view)

        self.previousMCA.setEnabled(False)
        self.nextMCA    .setEnabled(False)

        self.record = []
        self.record_i = 0

    def set_ROI_widgets(self):
        """
        Appends all ROI related fields to their related lists in Display for reference
        """
        for i in range(1, self.num_ROI + 1):
            self.ROI.append(self.findChild(QtWidgets.QCheckBox, f"ROI{i}"))
            self.start.append(self.findChild(QtWidgets.QLineEdit, f"start{i}"))
            self.end.append(self.findChild(QtWidgets.QLineEdit, f"end{i}"))
            self.counts.append(
                self.findChild(
                    QtWidgets.QLineEdit,
                    f"counts{i}"))
            self.lines.append(self.findChild(QtWidgets.QLineEdit, f"lines{i}"))
        return

    def ui_filename(self):
        return 'SSRL_MCA.ui'

    def ui_filepath(self):
        return path.join(path.dirname(path.realpath(__file__)),
                         self.ui_filename())

    def mouse_moved(self, point):
        """
        Tracks the mouse movement in the view
        """
        if (not self.waveform.sceneBoundingRect().contains(point)):
            return

        point_v = self.waveform.getViewBox().mapSceneToView(point)

        emin = int(point_v.x()) - 200
        emax = int(point_v.x()) + 200
        line_e = [
            ei for ei in self.energy if (
                (ei[0] > emin) and (
                    ei[0] < emax))]
        line_p = sorted(line_e, key=itemgetter(2))

        l_text = ""

        for ip in range(min(6, len(line_p))):
            if (ip > 0):
                l_text = l_text + ", "

            l_text = l_text + line_p[ip][1] + "-" + line_p[ip][2] + ": " +     \
                str(int(line_p[ip][0]))

        self.mouse_e.setText(str(int(point_v.x())))
        self.mouse_c.setText(str(int(point_v.y())))
        self.mouse_p.setText(l_text)

    def parse_args(self, args):
        parser = argparse.ArgumentParser()
        parser.add_argument(
            '--f',
            dest='filename',
            help='Input filename as string to be opened.')
        parser.add_argument(
            '--d',
            dest='data',
            help='Indicator that data should be read through stream.')
        parsed_args, _unknown_args = parser.parse_known_args(args)
        return parsed_args

    def full_view(self, *args, **kwargs):
        #       self.waveform.enableAutoRange()

        self.waveform.resetAutoRangeX()
        self.waveform.resetAutoRangeY()

    def change_source(self, *args, **kwargs):
        if (args[0] == 0):
            self.show_exposure()
        else:
            self.show_mca()
        #self.previousMCA.setEnabled(False)
        #self.nextMCA    .setEnabled(False)

    def show_exposure(self):
        """
        self.recordNum_l    .hide()
        self.recordNum      .hide()
        self.openFile       .hide()
        self.previousMCA    .hide()
        self.nextMCA        .hide()
        """

        self.exposure_l     .setEnabled(True)
        self.exposure       .setEnabled(True)
        self.exposureCount_l.show()
        self.exposureCount  .show()
        self.start_b        .show()
        self.stop_b         .show()
        return

    def show_mca(self):
        self.recordNum_l    .show()
        self.recordNum      .show()
        self.openFile       .show()
        self.previousMCA    .show()
        self.nextMCA        .show()
        self.start_b        .show()
        self.stop_b         .show()

        """
        self.exposure_l     .setEnabled(False)
        self.exposure       .setEnabled(False)
        self.exposureCount_l.hide()
        self.exposureCount  .hide()
        self.start_b        .hide()
        self.stop_b         .hide()
        """
        return

    def connect_data(self):
        self.epics.connect()
        return 

    def disconnect_data(self):
        if self.epics:
            self.epics.disconnect()
        return

    def live_data(self, new_waveform):
        self.record = new_waveform
        # TODO: this is going to be overwriting the entire list with only one
        # waveform? Shouldn't this just be appended so that we can go backwords and 
        #forwards as wished?

        self.handle_mca()

    def open_file(self, *args, **kwargs):
        """
        Opens the data file specified by the user and sends for MCA processing
        """
        fname = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open file", "", "Data files (*.dat);;All files (*.*)")
        # No file selected
        if (fname[0] == ""):
            return

        self.openFilename.setText(fname[0])
        # TODO: parses all at once and throws into one list of line by line, I imagine I'll need to just populate a line at a time as it comes in?
        with open(fname[0]) as f:
            self.record = [line.rstrip() for line in f]
        self.record_i = 0

        if (len(self.record) > 0):
            if (len(self.record) > 1):
                self.nextMCA.setEnabled(True)
            self.handle_mca()
        self.previousMCA.setEnabled(False)

    def previous_mca(self, *args, **kwargs):
        self.logger.debug("\nPrevious MCA ...{}".format(self.separator))
        self.record_i = self.record_i - 1
        if (self.record_i == 0):
            self.previousMCA.setEnabled(False)

        self.nextMCA.setEnabled(True)

        self.handle_mca()

    def next_mca(self, *args, **kwargs):
        self.logger.debug("\nNext MCA ...{}".format(self.separator))
        self.record_i = self.record_i + 1
        if (self.record_i == len(self.record) - 1):
            self.nextMCA.setEnabled(False)

        self.previousMCA.setEnabled(True)

        self.handle_mca()

    def find_peak(self, y_array):
        start = math.floor(int(self.start0.text()) / 10.)

        ret_i = []
        work_d = []
        for ri in range(self.num_ROI):
            if not (self.ROI[ri].isChecked()):
                continue
            try:
                xl = math.floor(int(self.start[ri].text()) / 10.)
                xr = math.floor(int(self.end[ri].text()) / 10.)
                points = xr - xl + 1
            except BaseException:
                continue

            self.logger.debug(
                "\nROI {0}\nxl - {1}\nxr - {2}{3}".format(ri, xl, xr, self.separator))
            if (points < 12):
                continue

            xl = xl - start
            xr = xr - start

            ypeak = max(y_array[xl:xr + 1])
            xpeak = y_array[xl:xr + 1].index(ypeak) + xl

            self.logger.debug(
                "\nFit 0: \nri - {0}\nxl - {1}, xr - {2}\nxpeak - {3}, ypeak - {4}{5}".format(
                    ri, xl, xr, xpeak, ypeak, self.separator))
            try:
                if (self.fitc == "Cauchy"):
                    fit, tmp = curve_fit(cauchy, list(range(
                        xl, xr + 1)), y_array[xl:xr + 1], p0=[ypeak, (xr + xl) / 2., (xr - xl) / 4.])
                else:
                    fit, tmp = curve_fit(gaussian, list(range(
                        xl, xr + 1)), y_array[xl:xr + 1], p0=[ypeak, (xr + xl) / 2., (xr - xl) / 4.])
                fit = list(fit)
            except BaseException:
                fit = []

            # try to fit 2 curves
            if (fit != []) and (((fit[1] - xl) / (xr - xl) < 0.35)
                                or ((fit[1] - xl) / (xr - xl) > 0.65)):
                try:
                    fit2, tmp = curve_fit(gaussian2, list(range(xl, xr + 1)), y_array[xl:xr + 1], p0=[
                                          fit[0], fit[1], fit[2], self.num_ROI, xl + xr - fit[1], (xr - xl) / 4.])
                    self.logger.debug(
                        "\nFit2: {}{}".format(
                            fit2, self.separator))
                    fit = fit2
                except BaseException:
                    pass

            self.logger.debug(
                "\nFit i: \nxl - {}, xr - {}\nfit - {}{}".format(xl, xr, fit, self.separator))

            ret_i.append([xl, xr, ypeak, xpeak, fit])

            work_d.append([xl, xr])

        work_i = sorted(work_d, key=itemgetter(0))

        work_l = []
        xi = 0
        for wi in range(len(work_i)):
            xl = work_i[wi][0]
            xr = work_i[wi][1]

            if (xl - xi >= 12):
                ymax = max(y_array[xi:xl])
                if (ymax >= 80):
                    work_l.append([ymax, xi, xl - 1])

            xi = xr + 1

        if (len(y_array) - xi >= 12):
            ymax = max(y_array[xi:])
            if (ymax >= 80):
                work_l.append([ymax, xi, len(y_array) - 1])

        ret_l = []
        while(len(work_l) > 0):
            work_i = sorted(work_l, key=itemgetter(0), reverse=True)
            work_l = work_i[1:]
            self.logger.debug(
                "\nWork List: {}{}".format(
                    work_i, self.separator))

            if (work_i[0][0] < 80):
                continue                  # counts too low

            ypeak = work_i[0][0]
            xmin = work_i[0][1]
            xmax = work_i[0][2]

            y = y_array[xmin:xmax + 1]
            xpeak = y.index(ypeak)

            # monotonically going up or down
            if ((ypeak == y[0]) and (min(y) == y[-1])) or                   \
               ((ypeak == y[-1]) and (min(y) == y[0])):
                continue

            # ending up
            xr = len(y) - 3
            while(xr >= self.num_ROI):
                slope, intercept = np.polyfit(range(3), y[xr:xr + 3], 1)
                if (slope < 0.):
                    break

                xr = xr - 3

            xr = xr + 3

            if (xr < 12):
                continue                  # less than 12 points left

            xl = 0                                               # starting down
            while(xl <= xr - 12):                # xr is the length of the new y
                slope, intercept = np.polyfit(range(3), y[xl:xl + 3], 1)
                if (slope > 0.):
                    break

                xl = xl + 3

            if (xr - xl < 12):
                continue               # less than 12 points left

            xmax = xmin + xr - 1
            xmin = xmin + xl

            y = y_array[xmin:xmax + 1]
            ypeak = max(y)
            xpeak = y.index(ypeak)

            if (ypeak < 80):
                continue                         # counts too low

            dx = 10
            smax = 0
            xl = xpeak - dx + 1
            while((xl >= 0) and (xl > xpeak - 100)):
                slope, intercept = np.polyfit(range(dx), y[xl:xl + dx], 1)

                if (slope > smax):
                    smax = slope
                elif (slope < 0.5):
                    if (dx == 3):
                        if ((y[xl] + y[xl + dx - 1]) < ypeak * 0.8) or  \
                           (((y[xl] + y[xl + dx - 1]) < ypeak) and (slope < 0)):
                            xl = xl + 3
                            break
                    elif (dx == 5):
                        dx = 3
                        xl = xl + 8
                        continue
                    else:
                        dx = 5
                        xl = xl + 15
                        continue

                if (xl >= dx):
                    xl = xl - dx
                else:
                    break

            dx = 10
            smin = 0
            xr = xpeak
            while((xr <= len(y) - dx) and (xr < xpeak + 100)):
                slope, intercept = np.polyfit(range(dx), y[xr:xr + dx], 1)

                if (slope < smin):
                    smin = slope
                elif (slope > -0.5):
                    if (dx == 3):
                        if ((y[xr] + y[xr + dx - 1]) < ypeak * 0.8) or  \
                           (((y[xr] + y[xr + dx - 1]) < ypeak) and (slope > 0)):
                            xr = xr - 3
                            break
                    elif (dx == 5):
                        dx = 3
                        xr = xr - 3
                        continue
                    else:
                        dx = 5
                        xr = xr - 5
                        continue

                xr = xr + dx

            try:
                if (self.fitc == "Cauchy"):
                    fit, tmp = curve_fit(cauchy, list(
                        range(xl, xr + 1)), y[xl:xr + 1], p0=[ypeak, (xr + xl) / 2., (xr - xl) / 4.])
                else:
                    fit, tmp = curve_fit(gaussian, list(
                        range(xl, xr + 1)), y[xl:xr + 1], p0=[ypeak, (xr + xl) / 2., (xr - xl) / 4.])
                fit = list(fit)
            except BaseException:
                fit = []

            ret_l.append([xl + xmin, xr + xmin, ypeak, xpeak + xmin, fit])

            self.logger.debug(
                "\nFit: \nxl+xmin - {}\nxr+xmin - {}\nfit - {}{}".format(
                    xl + xmin, xr + xmin, fit, self.separator))

            if (len(ret_i) + len(ret_l) == 10):
                break

            if (xl >= 12):
                work_l.append([max(y[:xl]), xmin, xmin + xl - 1])

            if (len(y) - xr >= 13):
                work_l.append([max(y[xr + 1:]), xmin + xr + 1, xmax])

        return ret_i + sorted(ret_l, key=itemgetter(0))

    def handle_mca(self):
        # TODO: 
        # Seems like it should have a different comparison
        if (self.dataSource.currentText() == "Live EPICS"):
            # Need to set self.record based on data coming in
            print("in live epics")
            items = self.record
            print(self.record)
        else:
            items = list(map(int, self.record[self.record_i].split()))
            print("other items")
            print(items)
        if not items.any():
            return

        start = math.floor(int(self.start0.text()) / 10.)
        end = math.ceil(int(self.end0  .text()) / 10.) + 1

        order = 3
        cutoff = 0.1
        B, A = signal.butter(order, cutoff, output='ba')

        y = signal.filtfilt(B, A, items[1:][start:end])
        y_array = list(y)

#       y_array = items[1:][start:end]

        ymax = max(y_array)
        sum0 = sum(y_array)

        self.curve.receiveXWaveform(
            np.array(
                list(
                    range(
                        start *
                        10,
                        end *
                        10,
                        10))))
        self.curve.receiveYWaveform(np.array(y_array))

        self.counts0.setText('{:.0f}'.format(sum0))

        ret_l = self.find_peak(y_array)
        self.logger.debug("\nret_l - {0}{1}".format(ret_l, self.separator))

        for il in range(self.num_ROI):
            if (len(ret_l) > il):
                self.croi[il].show()
            else:
                self.croi[il].hide()
                self.line[il * 2].hide()
                self.line[il * 2 + 1].hide()

                self.counts[il].setText("")
                self.lines[il].setText("")

                if not (self.ROI[il].isChecked()):
                    self.start[il].setText("")
                    self.end[il].setText("")

                continue

            start_i = start + ret_l[il][0]
            end_i = start + ret_l[il][1]

            if (start_i < 0):
                start_i = 0
            if (end_i > 2047):
                end_i = 2047

            end_i = end_i + 1

            y_array = items[1:][start_i:end_i]
            ysum = sum(y_array)

            self.counts[il].setText('{:.0f}'.format(ysum))

            self.croi[il].receiveXWaveform(
                np.array(
                    list(
                        range(
                            start_i *
                            10,
                            end_i *
                            10,
                            10))))

            l_text = ""
            fit = ret_l[il][4]
            if not(self.ROI[il].isChecked()):
                self.start[il].setText(str(10 * (start_i)))
                self.end[il].setText(str(10 * (end_i)))

                self.croi[il].receiveYWaveform(np.array(y_array))

                if (len(fit) > 2):
                    efit = 10 * (start + fit[1])
                    efit = 10 * (start + ret_l[il][3])
                    emin = efit - 10 * fit[2]
                    emax = efit + 10 * fit[2]
                else:
                    efit = 10 * (start + ret_l[il][3])
                    emin = efit * 0.99
                    emax = efit * 1.01

                if (ret_l[il][2] > 0.4 * ymax):
                    scale = 1.1
                else:
                    scale = 0.5

                self.line[il * 2].receiveXWaveform(np.array([efit, efit]))
                self.line[il * 2].receiveYWaveform(np.array([0, scale * ymax]))
                self.line[il * 2].show()
                self.line[il * 2 + 1].hide()

                line_e = [
                    ei for ei in self.energy if (
                        (ei[0] > emin) and (
                            ei[0] < emax))]
                line_p = sorted(line_e, key=itemgetter(2))
                self.logger.debug(
                    "\nEfit - {}\nline - {}{}".format(efit, line_p, self.separator))

                for ip in range(min(6, len(line_p))):
                    if (ip > 0):
                        l_text = l_text + ", "

                    l_text = l_text + line_p[ip][1] + "-" + line_p[ip][2]
            elif (len(fit) < 3):
                self.croi[il].receiveYWaveform(np.array(y_array))

                efit = 10 * (start + ret_l[il][3])

                if (ret_l[il][2] > 0.4 * ymax):
                    scale = 1.1
                else:
                    scale = 0.5

                self.line[il * 2].receiveXWaveform(np.array([efit, efit]))
                self.line[il * 2].receiveYWaveform(np.array([0, scale * ymax]))
                self.line[il * 2].show()
                self.line[il * 2 + 1].hide()

                l_text = "Failed to fit"
            elif (len(fit) == 6):
                self.croi[il].receiveYWaveform(
                    gaussian2(
                        list(
                            range(
                                start_i * 10,
                                end_i * 10,
                                10)),
                        fit[0],
                        (fit[1] + start) * 10,
                        fit[2] * 10,
                        fit[3],
                        (fit[4] + start) * 10,
                        fit[5] * 10))

                if (fit[0] >= 50):
                    efit = 10 * (start + fit[1])

                    if (fit[0] > 0.4 * ymax):
                        scale = 1.1
                    else:
                        scale = 0.5

                    self.line[il * 2].receiveXWaveform(np.array([efit, efit]))
                    self.line[il *
                              2].receiveYWaveform(np.array([0, scale * ymax]))
                    self.line[il * 2].show()
                else:
                    self.line[il * 2].hide()

                if (fit[3] >= 50):
                    efit = 10 * (start + fit[4])

                    if (fit[3] > 0.4 * ymax):
                        scale = 1.1
                    else:
                        scale = 0.5

                    self.line[il * 2 +
                              1].receiveXWaveform(np.array([efit, efit]))
                    self.line[il * 2 +
                              1].receiveYWaveform(np.array([0, scale * ymax]))
                    self.line[il * 2 + 1].show()
                else:
                    self.line[il * 2 + 1].hide()

                l_text = str(int(fit[0])) + " " + str(int((start + fit[1]) * 10)) + " " + str(int(fit[2] * 10)) + "; " + \
                    str(int(fit[3])) + " " + str(int((start + fit[4])
                                                     * 10)) + " " + str(int(fit[5] * 10))
            elif (self.fitc == "Cauchy"):
                self.croi[il].receiveYWaveform(
                    cauchy(
                        list(
                            range(
                                start_i * 10,
                                end_i * 10,
                                10)),
                        fit[0] * 10,
                        (fit[1] + start) * 10,
                        fit[2] * 10))

                efit = 10 * (start + fit[1])

                if (ret_l[il][2] > 0.4 * ymax):
                    scale = 1.1
                else:
                    scale = 0.5

                self.line[il * 2].receiveXWaveform(np.array([efit, efit]))
                self.line[il * 2].receiveYWaveform(np.array([0, scale * ymax]))
                self.line[il * 2].show()
                self.line[il * 2 + 1].hide()

                l_text = str(
                    int(fit[0])) + " " + str(int((start + fit[1]) * 10)) + " " + str(int(fit[2] * 10))
            else:
                self.croi[il].receiveYWaveform(
                    gaussian(
                        list(
                            range(
                                start_i * 10,
                                end_i * 10,
                                10)),
                        fit[0],
                        (fit[1] + start) * 10,
                        fit[2] * 10))

                efit = 10 * (start + fit[1])

                if (ret_l[il][2] > 0.4 * ymax):
                    scale = 1.1
                else:
                    scale = 0.5

                self.line[il * 2].receiveXWaveform(np.array([efit, efit]))
                self.line[il * 2].receiveYWaveform(np.array([0, scale * ymax]))
                self.line[il * 2].show()
                self.line[il * 2 + 1].hide()

                l_text = str(
                    int(fit[0])) + " " + str(int((start + fit[1]) * 10)) + " " + str(int(fit[2] * 10))

            self.lines[il].setText(l_text)

        self.recordNum.setText(str(items[0]))
