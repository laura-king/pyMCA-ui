import math, time
import numpy                   as np
import scipy.signal            as signal

from qtpy.QtCore           import Qt
from qtpy                  import QtGui
from qtpy.QtWidgets        import QFileDialog
from pydm                  import Display
from pydm.widgets.channel  import PyDMChannel
from scipy.optimize        import curve_fit
from operator              import itemgetter
from os                    import path, popen
import pprint

element  = {}
emission = {}
energy   = []

def build_dic():
    global energy

    line = popen( "cat XRay_Emission_Lines.txt" ).readlines()

    for il in range( len(line) ):
        if ( len(line[il]) <= 1 ): continue

        xxx = line[il].split()

        try:
            idx = int(xxx[0])
        except:
            continue

        symbol = xxx[1] + xxx[0]

        try:
            ka1 = float(xxx[2])
        except:
            ka1 = 0.

        try:
            ka2 = float(xxx[3])
        except:
            ka2 = 0.

        try:
            kb1 = float(xxx[4])
        except:
            kb1 = 0.

        try:
            la1 = float(xxx[5])
        except:
            la1 = 0.

        try:
            la2 = float(xxx[6])
        except:
            la2 = 0.

        try:
            lb1 = float(xxx[7])
        except:
            lb1 = 0.

        try:
            lb2 = float(xxx[8])
        except:
            lb2 = 0.

        try:
            lg1 = float(xxx[9])
        except:
            lg1 = 0.

        try:
            ma1 = float(xxx[10])
        except:
            ma1 = 0.

        element_d = {}

        if ( ka1 != 0 ) and ( ka2 != 0 ):
            if ( (ka1 - ka2) > -30 ) and ( (ka1 - ka2) < 30 ):
                element_d[        "Ka" ] = (ka1 + ka2) / 2.

                emission [symbol+"-Ka" ] = (ka1 + ka2) / 2.
            else:
                element_d[        "Ka1"] =  ka1
                element_d[        "Ka2"] =  ka2

                emission [symbol+"-Ka1"] =  ka1
                emission [symbol+"-Ka2"] =  ka2
        else:
            if ( ka1 != 0 ):
                element_d[        "Ka1"] =  ka1

                emission [symbol+"-Ka1"] =  ka1

            if ( ka2 != 0 ):
                element_d[        "Ka2"] =  ka2

                emission [symbol+"-Ka2"] =  ka2

        if ( kb1 != 0 ):
            element_d[        "Kb1"] =  kb1

            emission [symbol+"-Kb1"] =  kb1

        if ( la1 != 0 ) and ( la2 != 0 ):
            if ( (la1 - la2) > -30 ) and ( (la1 - la2) < 30 ):
                element_d[        "La" ] = (la1 + la2) / 2.

                emission [symbol+"-La" ] = (la1 + la2) / 2.
            else:
                element_d[        "La1"] =  la1
                element_d[        "La2"] =  la2

                emission [symbol+"-La1"] =  la1
                emission [symbol+"-La2"] =  la2
        else:
            if ( la1 != 0 ):
                element_d[        "La1"] =  la1

                emission [symbol+"-La1"] =  la1

            if ( la2 != 0 ):
                element_d[        "La2"] =  la2

                emission [symbol+"-La2"] =  la2

        if ( lb1 != 0 ) and ( lb2 != 0 ):
            if ( (lb1 - lb2) > -30 ) and ( (lb1 - lb2) < 30 ):
                element_d[        "Lb" ] = (lb1 + lb2) / 2.

                emission [symbol+"-Lb" ] = (lb1 + lb2) / 2.
            else:
                element_d[        "Lb1"] =  lb1
                element_d[        "Lb2"] =  lb2

                emission [symbol+"-Lb1"] =  lb1
                emission [symbol+"-Lb2"] =  lb2
        else:
            if ( lb1 != 0 ):
                element_d[        "Lb1"] =  lb1

                emission [symbol+"-Lb1"] =  lb1

            if ( lb2 != 0 ):
                element_d[        "Lb2"] =  lb2

                emission [symbol+"-Lb2"] =  lb2

        if ( lg1 != 0 ):
            element_d[        "Lg1"] =  lg1

            emission [symbol+"-Lg1"] =  lg1

        if ( ma1 != 0 ):
            element_d[        "Ma1"] =  ma1

            emission [symbol+"-Ma1"] =  ma1

        element[symbol] = element_d
    pprint.pprint(emission)
    pprint.pprint(element_d)

    energy_i = []
    for key in emission.keys():
        xxx = key.split( "-" )
        energy_i.append( [emission[key], xxx[0], xxx[1]] )

    energy = sorted( energy_i, key=itemgetter(0) )

def gaussian (x, amplitude, mean, sigma):
    return amplitude * np.exp(-((x-mean)/sigma)**2/2.)

def gaussian2(x, ampl1, mean1, sigma1, ampl2, mean2, sigma2):
    return ampl1 * np.exp(-((x-mean1)/sigma1)**2/2.) +                         \
           ampl2 * np.exp(-((x-mean2)/sigma2)**2/2.)

def cauchy(x, amplitude, mean, gamma):
    return amplitude * gamma / ((x-mean)*(x-mean) + gamma*gamma)

class MCADisplay( Display ):
    def __init__(self, parent=None, args=None, macros=None):
        super(MCADisplay, self).__init__(parent=parent, args=args, macros=macros)

        build_dic()

#       self.waveform.enableCrosshair(True, 0, 0)
#       self.waveform.crosshair_position_updated.connect    (self.cross_info)

        self.waveform.plotItem.scene().sigMouseMoved.connect(self.mouseMoved)

        self.waveform.setXLabels(["Energy (eV)"])
        self.waveform.setYLabels(["Count"      ])

        self.waveform.addChannel(None, None, name="Full",   color="white")

        self.waveform.addChannel(None, None, name="ROI1",   color="red"  ,     \
                                             lineWidth=2)
        self.waveform.addChannel(None, None, name="ROI2",   color="green",     \
                                             lineWidth=2)
        self.waveform.addChannel(None, None, name="ROI3",   color="blue" ,     \
                                             lineWidth=2)
        self.waveform.addChannel(None, None, name="ROI4",   color="red"  ,     \
                                             lineWidth=2)
        self.waveform.addChannel(None, None, name="ROI5",   color="green",     \
                                             lineWidth=2)
        self.waveform.addChannel(None, None, name="ROI6",   color="blue" ,     \
                                             lineWidth=2)
        self.waveform.addChannel(None, None, name="ROI7",   color="red"  ,     \
                                             lineWidth=2)
        self.waveform.addChannel(None, None, name="ROI8",   color="green",     \
                                             lineWidth=2)
        self.waveform.addChannel(None, None, name="ROI9",   color="blue" ,     \
                                             lineWidth=2)

        self.waveform.addChannel(None, None, name="Line01", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line02", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line03", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line04", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line05", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line06", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line07", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line08", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line09", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line11", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line12", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line13", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line14", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line15", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line16", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line17", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line18", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)
        self.waveform.addChannel(None, None, name="Line19", color="white",     \
                                             lineWidth=2, lineStyle=Qt.DashLine)

        self.curve = self.waveform._curves[ 0];

        self.croi  = []
        self.croi.append( self.waveform._curves[ 1] );
        self.croi.append( self.waveform._curves[ 2] );
        self.croi.append( self.waveform._curves[ 3] );
        self.croi.append( self.waveform._curves[ 4] );
        self.croi.append( self.waveform._curves[ 5] );
        self.croi.append( self.waveform._curves[ 6] );
        self.croi.append( self.waveform._curves[ 7] );
        self.croi.append( self.waveform._curves[ 8] );
        self.croi.append( self.waveform._curves[ 9] );

        self.line  = []
        self.line.append( self.waveform._curves[10] );
        self.line.append( self.waveform._curves[11] );
        self.line.append( self.waveform._curves[12] );
        self.line.append( self.waveform._curves[13] );
        self.line.append( self.waveform._curves[14] );
        self.line.append( self.waveform._curves[15] );
        self.line.append( self.waveform._curves[16] );
        self.line.append( self.waveform._curves[17] );
        self.line.append( self.waveform._curves[18] );
        self.line.append( self.waveform._curves[19] );
        self.line.append( self.waveform._curves[20] );
        self.line.append( self.waveform._curves[21] );
        self.line.append( self.waveform._curves[22] );
        self.line.append( self.waveform._curves[23] );
        self.line.append( self.waveform._curves[24] );
        self.line.append( self.waveform._curves[25] );
        self.line.append( self.waveform._curves[26] );
        self.line.append( self.waveform._curves[27] );

        self.ROI    = []
        self.ROI   .append( self.ROI1   );
        self.ROI   .append( self.ROI2   );
        self.ROI   .append( self.ROI3   );
        self.ROI   .append( self.ROI4   );
        self.ROI   .append( self.ROI5   );
        self.ROI   .append( self.ROI6   );
        self.ROI   .append( self.ROI7   );
        self.ROI   .append( self.ROI8   );
        self.ROI   .append( self.ROI9   );

        self.start  = []
        self.start .append( self.start1 );
        self.start .append( self.start2 );
        self.start .append( self.start3 );
        self.start .append( self.start4 );
        self.start .append( self.start5 );
        self.start .append( self.start6 );
        self.start .append( self.start7 );
        self.start .append( self.start8 );
        self.start .append( self.start9 );

        self.end    = []
        self.end   .append( self.end1   );
        self.end   .append( self.end2   );
        self.end   .append( self.end3   );
        self.end   .append( self.end4   );
        self.end   .append( self.end5   );
        self.end   .append( self.end6   );
        self.end   .append( self.end7   );
        self.end   .append( self.end8   );
        self.end   .append( self.end9   );

        self.counts = []
        self.counts.append( self.counts1 );
        self.counts.append( self.counts2 );
        self.counts.append( self.counts3 );
        self.counts.append( self.counts4 );
        self.counts.append( self.counts5 );
        self.counts.append( self.counts6 );
        self.counts.append( self.counts7 );
        self.counts.append( self.counts8 );
        self.counts.append( self.counts9 );

        self.lines  = []
        self.lines .append( self.lines1  );
        self.lines .append( self.lines2  );
        self.lines .append( self.lines3  );
        self.lines .append( self.lines4  );
        self.lines .append( self.lines5  );
        self.lines .append( self.lines6  );
        self.lines .append( self.lines7  );
        self.lines .append( self.lines8  );
        self.lines .append( self.lines9  );

        if ( macros != None ) and ( "FIT" in macros ):
            if ( macros["FIT"].lower() == "cauchy" ): self.fitc = "Cauchy"
            else:                                     self.fitc = "Gaussian"
        else:
            self.fitc = "Gaussian"

        if ( macros != None ) and ( "DEVICE" in macros ):
            self.dataSource.addItem( "Live EPICS" )

            epics = PyDMChannel( address="ca://"+
                                 macros["DEVICE"]+":ARR1:ArrayData",
                                 value_slot=self.live_data )
            epics.connect()

            self.recordNum_l    .hide()
            self.recordNum      .hide()
            self.openFile       .hide()
            self.previousMCA    .hide()
            self.nextMCA        .hide()

            self.exposure_l     .setEnabled( True  )
            self.exposure       .setEnabled( True  )
            self.exposureCount_l.show()
            self.exposureCount  .show()
            self.start_b        .show()
            self.stop_b         .show()
        else:
            self.recordNum_l    .show()
            self.recordNum      .show()
            self.openFile       .show()
            self.previousMCA    .show()
            self.nextMCA        .show()

            self.exposure_l     .setEnabled( False )
            self.exposure       .setEnabled( False )
            self.exposureCount_l.hide()
            self.exposureCount  .hide()
            self.start_b        .hide()
            self.stop_b         .hide()

        self.dataSource.addItem        ( "Playback" )
        self.dataSource.setCurrentIndex( 0 )
        self.dataSource.currentIndexChanged.connect( self.change_source )

        self.openFile   .clicked           .connect( self.open_file     )
        self.previousMCA.clicked           .connect( self.previous_mca  )
        self.nextMCA    .clicked           .connect( self.next_mca      )
        self.fullView   .clicked           .connect( self.full_view     )

        self.previousMCA.setEnabled( False )
        self.nextMCA    .setEnabled( False )

        self.record   = []
        self.record_i = 0

    def ui_filename(self):
        return 'SSRL_MCA.ui'

    def ui_filepath(self):
        return path.join(path.dirname(path.realpath(__file__)),                \
                         self.ui_filename())

    def mouseMoved(self, point):
        if ( not self.waveform.sceneBoundingRect().contains(point) ):
            return

        point_v = self.waveform.getViewBox().mapSceneToView(point)

        emin    = int(point_v.x()) - 200
        emax    = int(point_v.x()) + 200
        line_e  = [ ei for ei in energy if ((ei[0] > emin) and (ei[0] < emax)) ]
        line_p  = sorted( line_e, key=itemgetter(2) )

        l_text  = ""
        for ip in range( min(6, len(line_p)) ):
            if ( ip > 0 ): l_text = l_text + ", "

            l_text = l_text + line_p[ip][1] + "-" + line_p[ip][2] + ": " +     \
                     str(int(line_p[ip][0]))

        self.mouse_e.setText( str(int(point_v.x())) )
        self.mouse_c.setText( str(int(point_v.y())) )
        self.mouse_p.setText( l_text                )

    def full_view(self, *args, **kwargs):
#       self.waveform.enableAutoRange()

        self.waveform.resetAutoRangeX()
        self.waveform.resetAutoRangeY()

    def change_source(self, *args, **kwargs):
        if ( args[0] == 0 ):
            self.recordNum_l    .hide()
            self.recordNum      .hide()
            self.openFile       .hide()
            self.previousMCA    .hide()
            self.nextMCA        .hide()

            self.exposure_l     .setEnabled( True  )
            self.exposure       .setEnabled( True  )
            self.exposureCount_l.show()
            self.exposureCount  .show()
            self.start_b        .show()
            self.stop_b         .show()
        else:
            self.recordNum_l    .show()
            self.recordNum      .show()
            self.openFile       .show()
            self.previousMCA    .show()
            self.nextMCA        .show()

            self.exposure_l     .setEnabled( False )
            self.exposure       .setEnabled( False )
            self.exposureCount_l.hide()
            self.exposureCount  .hide()
            self.start_b        .hide()
            self.stop_b         .hide()

        self.previousMCA.setEnabled( False )
        self.nextMCA    .setEnabled( False )

    def live_data(self, new_waveform):
        self.record = new_waveform

        self.handle_mca()

    def open_file(self, *args, **kwargs):
        fname = QFileDialog.getOpenFileName(self, "Open file", "", "Data files (*.dat);;All files (*.*)" )

        if ( fname[0] != "" ): self.record = popen("cat "+fname[0]).readlines()
        else:                  self.record = []

        self.record_i = 0
        if ( len(self.record) > 0 ):
            if ( len(self.record) > 1 ): self.nextMCA.setEnabled( True )

            self.handle_mca()

        self.previousMCA.setEnabled( False )

    def previous_mca(self, *args, **kwargs):
        print('Previous MCA ...')

        self.record_i = self.record_i - 1
        if ( self.record_i == 0 ): self.previousMCA.setEnabled( False )

        self.nextMCA.setEnabled( True )

        self.handle_mca()

    def next_mca(self, *args, **kwargs):
        print('Next MCA ...')

        self.record_i = self.record_i + 1
        if ( self.record_i == len(self.record)-1 ): self.nextMCA.setEnabled( False )

        self.previousMCA.setEnabled( True )

        self.handle_mca()

    def find_peak(self, y_array):
        start  = math.floor(int(self.start0.text()) / 10.)

        ret_i  = []
        work_d = []
        for ri in range(9):
            if ( self.ROI[ri].checkState() != 2 ): continue

            print( "ROI: ", ri )
            try:
                xl     = math.floor(int(self.start[ri].text()) / 10.)
                xr     = math.floor(int(self.end  [ri].text()) / 10.)
                points = xr - xl + 1
            except:
                continue

            print( "ROI: ", ri, xl, xr )
            if ( points < 12 ): continue

            xl    = xl - start
            xr    = xr - start

            ypeak = max(y_array[xl:xr+1])
            xpeak = y_array[xl:xr+1].index(ypeak) + xl

            print( "Fit0 ", ri, xl, xr, xpeak, ypeak )
            try:
                if ( self.fitc == "Cauchy" ):
                    fit, tmp = curve_fit(cauchy,    list(range(xl, xr+1)),     \
                                                    y_array[xl:xr+1],          \
                                                    p0=[ypeak, (xr+xl)/2., (xr-xl)/4.])
                else:
                    fit, tmp = curve_fit(gaussian,  list(range(xl, xr+1)),     \
                                                    y_array[xl:xr+1],          \
                                                    p0=[ypeak, (xr+xl)/2., (xr-xl)/4.])

#                   # try to fit 2 gaussians
#                   fit2,tmp = curve_fit(gaussian2, list(range(xl, xr+1)),     \
#                                                   y_array[xl:xr+1],          \
#                                                   p0=[fit[0], fit[1], fit[2], 9, xl+xr-fit[1], (xr-xl)/4.])
#                   print( "Fit2: ", fit2 )
#                   if ( len(fit2) > 0 ): fit = fit2

                fit = list(fit)
            except:
                fit = []

            # try to fit 2 curves
            if ( fit != [] ) and                                               \
               ( ((fit[1]-xl)/(xr-xl) < 0.35) or ((fit[1]-xl)/(xr-xl) > 0.65) ):
                try:
                    fit2, tmp = curve_fit(gaussian2, list(range(xl, xr+1)),    \
                                                     y_array[xl:xr+1],         \
                                                     p0=[fit[0], fit[1], fit[2], 9, xl+xr-fit[1], (xr-xl)/4.])
                    print( "Fit2: ", fit2 )
                    fit = fit2
                except:
                    pass

#           if ( len(fit) == 6 ) and ( fit[3] > fit[0] ):
#               fit = fit[3:] + fit[:3]

            print( "Fit i: ", xl, xr, fit )
            ret_i.append( [ xl, xr, ypeak, xpeak, fit ] )

            work_d.append([xl, xr])

        work_i = sorted( work_d, key=itemgetter(0) )

        work_l = []
        xi     = 0
        for wi in range( len(work_i) ):
            xl = work_i[wi][0]
            xr = work_i[wi][1]

            if ( xl-xi >= 12 ):
                ymax = max(y_array[xi:xl])
                if ( ymax >= 80 ): work_l.append([ymax, xi, xl-1])

            xi = xr + 1

        if ( len(y_array)-xi >= 12 ):
            ymax = max(y_array[xi:])
            if ( ymax >= 80 ): work_l.append([ymax, xi, len(y_array)-1])

        ret_l  = []
        while( len(work_l) > 0 ):
            work_i = sorted( work_l, key=itemgetter(0), reverse=True )
            work_l = work_i[1:]
            print( "Work List: ", work_i )

            if ( work_i[0][0] < 80 ): continue                  # counts too low

            ypeak  = work_i[0][0]
            xmin   = work_i[0][1]
            xmax   = work_i[0][2]

            y      = y_array[xmin:xmax+1]
            xpeak  = y.index(ypeak)

            # monotonically going up or down
            if ( (ypeak == y[ 0]) and (min(y) == y[-1]) ) or                   \
               ( (ypeak == y[-1]) and (min(y) == y[ 0]) ): continue

            xr = len(y) - 3                                          # ending up
            while( xr >= 9 ):
                slope, intercept = np.polyfit( range(3), y[xr:xr+3], 1 )
                if ( slope < 0. ): break

                xr = xr - 3

            xr = xr + 3

            if ( xr < 12 ): continue                  # less than 12 points left

            xl = 0                                               # starting down
            while( xl <= xr-12 ):                # xr is the length of the new y
                slope, intercept = np.polyfit( range(3), y[xl:xl+3], 1 )
                if ( slope > 0. ): break

                xl = xl + 3

            if ( xr-xl < 12 ): continue               # less than 12 points left

            xmax  = xmin + xr - 1
            xmin  = xmin + xl

            y     = y_array[xmin:xmax+1]
            ypeak = max(y)
            xpeak = y.index(ypeak)

            if ( ypeak < 80 ): continue                         # counts too low

            dx    = 10
            smax  =  0
            xl    = xpeak - dx + 1
            while( (xl >= 0) and (xl > xpeak - 100) ):
                slope, intercept = np.polyfit( range(dx), y[xl:xl+dx], 1 )

                if   ( slope > smax ): smax = slope
                elif ( slope < 0.5  ):
                    if   ( dx == 3 ):
                        if ( (y[xl]+y[xl+dx-1]) < ypeak*0.8              ) or  \
                           (((y[xl]+y[xl+dx-1]) < ypeak) and (slope < 0) ):
                            xl = xl +  3
                            break
                    elif ( dx == 5 ):
                        dx = 3
                        xl = xl +  8
                        continue
                    else:
                        dx = 5
                        xl = xl + 15
                        continue

                if ( xl >= dx ): xl = xl - dx
                else:            break

            dx    = 10
            smin  =  0
            xr    = xpeak
            while( (xr <= len(y)-dx) and (xr < xpeak + 100) ):
                slope, intercept = np.polyfit( range(dx), y[xr:xr+dx], 1 )

                if   ( slope < smin ): smin = slope
                elif ( slope > -0.5 ):
                    if   ( dx == 3 ):
                        if ( (y[xr]+y[xr+dx-1]) < ypeak*0.8              ) or  \
                           (((y[xr]+y[xr+dx-1]) < ypeak) and (slope > 0) ):
                            xr  = xr -  3
                            break
                    elif ( dx == 5 ):
                        dx = 3
                        xr = xr -  3
                        continue
                    else:
                        dx = 5
                        xr = xr -  5
                        continue

                xr = xr + dx

            try:
                if ( self.fitc == "Cauchy" ):
                    fit, tmp = curve_fit(cauchy,   list(range(xl, xr+1)),      \
                                                   y[xl:xr+1],                 \
                                                   p0=[ypeak, (xr+xl)/2., (xr-xl)/4.])
                else:
                    fit, tmp = curve_fit(gaussian, list(range(xl, xr+1)),      \
                                                   y[xl:xr+1],                 \
                                                   p0=[ypeak, (xr+xl)/2., (xr-xl)/4.])
                fit = list(fit)
            except:
                fit = []

            ret_l.append( [ xl+xmin, xr+xmin, ypeak, xpeak+xmin, fit ] )

            print( "Fit: ", xl+xmin, xr+xmin, fit )

            if ( len(ret_i)+len(ret_l) == 10 ): break

            if ( xl          >= 12 ):
                work_l.append( [max(y[:xl]  ), xmin,      xmin+xl-1] )

            if ( len(y) - xr >= 13 ):
                work_l.append( [max(y[xr+1:]), xmin+xr+1, xmax     ] )

        return ret_i + sorted( ret_l, key=itemgetter(0) )

    def handle_mca(self):
        if ( self.dataSource.currentText() == "Live EPICS" ):
            items = self.record
        else:
            items = list(map(int, self.record[self.record_i].split()))

        start   = math.floor(int(self.start0.text()) / 10.)
        end     = math.ceil (int(self.end0  .text()) / 10.) + 1
        mid     = (start + end) / 2

        order   = 3
        cutoff  = 0.1
        B, A    = signal.butter  ( order, cutoff, output='ba' )

        y       = signal.filtfilt( B, A, items[1:][start:end] )
        y_array = list( y )

#       y_array = items[1:][start:end]

        ymax    = max(y_array)
        sum0    = sum(y_array)

        self.curve.receiveXWaveform(np.array(list(range(start*10, end*10, 10))))
        self.curve.receiveYWaveform(np.array(y_array                          ))

        self.counts0.setText('{:.0f}'.format(sum0))

        ret_l = self.find_peak(y_array)
        print( ret_l )

        for il in range( 9 ):
            if ( len(ret_l) > il ):
                self.croi  [il    ].show()
            else:
                self.croi  [il    ].hide()
                self.line  [il*2  ].hide()
                self.line  [il*2+1].hide()

                self.counts[il    ].setText("")
                self.lines [il    ].setText("")

                if ( self.ROI[il].checkState() != 2 ):
                    self.start [il].setText("")
                    self.end   [il].setText("")

                continue

            start_i = start + ret_l[il][0]
            end_i   = start + ret_l[il][1]

            if (start_i < 0   ): start_i =    0
            if (end_i   > 2047): end_i   = 2047

            end_i   = end_i + 1

            y_array = items[1:][start_i:end_i]
            ysum    = sum(y_array)

            self.counts[il].setText('{:.0f}'.format(ysum))

            self.croi  [il].receiveXWaveform(np.array(list(range(start_i*10, end_i*10, 10))))

            l_text = ""
            fit    = ret_l[il][4]
            if ( self.ROI[il].checkState() != 2 ):
                self.start[il].setText(str(10*(start_i)))
                self.end  [il].setText(str(10*(end_i  )))

                self.croi [il].receiveYWaveform(np.array(y_array))

                if ( len(fit) > 2 ):
                    efit = 10 * ( start + fit[1] )
                    efit = 10 * ( start + ret_l[il][3] )
                    emin = efit - 10 * fit[2]
                    emax = efit + 10 * fit[2]
                else:
                    efit = 10 * ( start + ret_l[il][3] )
                    emin = efit * 0.99
                    emax = efit * 1.01

                if ( ret_l[il][2] > 0.4*ymax ): scale = 1.1
                else:                           scale = 0.5

                self.line[il*2  ].receiveXWaveform(np.array([efit, efit      ]))
                self.line[il*2  ].receiveYWaveform(np.array([0,    scale*ymax]))
                self.line[il*2  ].show()
                self.line[il*2+1].hide()

                line_e = [ ei for ei in energy if ((ei[0] > emin) and (ei[0] < emax)) ]
                line_p = sorted( line_e, key=itemgetter(2) )
                print( efit, line_p )

                for ip in range( min(6, len(line_p)) ):
                    if ( ip > 0 ): l_text = l_text + ", "

                    l_text = l_text + line_p[ip][1] + "-" + line_p[ip][2]
            elif ( len(fit) <  3 ):
                self.croi[il].receiveYWaveform(np.array(y_array))

                efit = 10 * ( start + ret_l[il][3] )

                if ( ret_l[il][2] > 0.4*ymax ): scale = 1.1
                else:                           scale = 0.5

                self.line[il*2  ].receiveXWaveform(np.array([efit, efit      ]))
                self.line[il*2  ].receiveYWaveform(np.array([0,    scale*ymax]))
                self.line[il*2  ].show()
                self.line[il*2+1].hide()

                l_text = "Failed to fit"
            elif ( len(fit) == 6 ):
                self.croi[il].receiveYWaveform(gaussian2(list(range(start_i*10, end_i*10, 10)), fit[0]   , (fit[1]+start)*10, fit[2]*10, fit[3]   , (fit[4]+start)*10, fit[5]*10))

                if ( fit[0] >= 50 ):
                    efit = 10 * ( start + fit[1] )

                    if ( fit[0] > 0.4*ymax ): scale = 1.1
                    else:                     scale = 0.5

                    self.line[il*2  ].receiveXWaveform(np.array([efit, efit      ]))
                    self.line[il*2  ].receiveYWaveform(np.array([0,    scale*ymax]))
                    self.line[il*2  ].show()
                else:
                    self.line[il*2  ].hide()

                if ( fit[3] >= 50 ):
                    efit = 10 * ( start + fit[4] )

                    if ( fit[3] > 0.4*ymax ): scale = 1.1
                    else:                     scale = 0.5

                    self.line[il*2+1].receiveXWaveform(np.array([efit, efit      ]))
                    self.line[il*2+1].receiveYWaveform(np.array([0,    scale*ymax]))
                    self.line[il*2+1].show()
                else:
                    self.line[il*2+1].hide()

                l_text = str(int(fit[0])) + " " + str(int((start+fit[1])*10)) + " " + str(int(fit[2]*10)) + "; " + \
                         str(int(fit[3])) + " " + str(int((start+fit[4])*10)) + " " + str(int(fit[5]*10))
            elif ( self.fitc == "Cauchy" ):
                self.croi[il].receiveYWaveform(cauchy   (list(range(start_i*10, end_i*10, 10)), fit[0]*10, (fit[1]+start)*10, fit[2]*10))

                efit = 10 * ( start + fit[1] )

                if ( ret_l[il][2] > 0.4*ymax ): scale = 1.1
                else:                           scale = 0.5

                self.line[il*2  ].receiveXWaveform(np.array([efit, efit      ]))
                self.line[il*2  ].receiveYWaveform(np.array([0,    scale*ymax]))
                self.line[il*2  ].show()
                self.line[il*2+1].hide()

                l_text = str(int(fit[0])) + " " + str(int((start+fit[1])*10)) + " " + str(int(fit[2]*10))
            else:
                self.croi[il].receiveYWaveform(gaussian (list(range(start_i*10, end_i*10, 10)), fit[0]   , (fit[1]+start)*10, fit[2]*10))

                efit = 10 * ( start + fit[1] )

                if ( ret_l[il][2] > 0.4*ymax ): scale = 1.1
                else:                           scale = 0.5

                self.line[il*2  ].receiveXWaveform(np.array([efit, efit      ]))
                self.line[il*2  ].receiveYWaveform(np.array([0,    scale*ymax]))
                self.line[il*2  ].show()
                self.line[il*2+1].hide()

                l_text = str(int(fit[0])) + " " + str(int((start+fit[1])*10)) + " " + str(int(fit[2]*10))

            self.lines[il].setText( l_text )

        self.recordNum.setText(str(items[0]))

