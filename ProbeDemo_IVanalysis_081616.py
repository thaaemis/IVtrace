from scipy import signal
import numpy as np
import os
from pylab import *
from IVclass import IVobj

def main():
    cwd = os.getcwd()
    
    path = "P:/departments/HTX/bkraus/ExpCampaign_EPionFlow_Aug2016/ProbeDemo_GodyakTest_081616/"
    # os.chdir(path)
    files = os.listdir(path)
    # os.chdir(cwd)
    txtFiles = [x for x in files if x[-4:] == '.txt']
    godyakFiles = [x for x in txtFiles if x[:7] != 'RawData']
    godyakFiles.remove('datalog.txt')
    IVlist = []
    for x in godyakFiles:
        ind = x[x.find('W')+1:]
        ind = ind[:ind.find('.txt')]
        ind = int(ind)
        thisIV = IVobj(x,path,flow = 0., B = 0, P = 0, V = 0,
                            I = 0.0, pos=0, ind=ind)
        # subplot(2,1,1)

        #     if ((thisIV.filter == 2) or (thisIV.filter == 2)) and (thisIV.sweepRate == 2):
        figure(1)
        isLabelonPlot = unicode(thisIV.label) in gca().get_legend_handles_labels()[1]
        thisIV.plotTrace(raw=False, deriv1 = False, deriv2 = False, label=not(isLabelonPlot))
        thisIV.plotTrace(raw=True, deriv1 = False, deriv2 = False, label=not(isLabelonPlot))
        
        figure(2)
        thisIV.plotEEDF()
        # print thisIV.ind, thisIV.B, thisIV.filter, thisIV.sweepRate, thisIV.Vf, thisIV.Vp, thisIV.ionSat
        thisIV.findTe(2,5)

    figure(1)
    manualVmeas = [0, 1, 2, 3, 4, 5, -5, 10]
    manualImeas = [-.0145, -.0074, .0034, .0221, .0552, .1139, -.0343, .886]
    manualV = [manualVmeas[i] - manualImeas[i] for i in range(len(manualVmeas))]
    manualI = [x/108. for x in manualImeas]
    plot(manualV, manualI, 'k*', markersize=20, label="Manual Points")
    legend(loc="best")
    text(6.5, 1e-4, 'Probe Demo Test',fontsize=20)
    
    show()
   

main()
