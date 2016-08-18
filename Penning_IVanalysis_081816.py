from scipy import signal
import numpy as np
import os
from pylab import *
from IVclass import IVobj

def main():
    cwd = os.getcwd()
    
    path = "P:/departments/HTX/bkraus/ExpCampaign_EPionFlow_Aug2016/LPEP_PenningProbeCheck_081616/"
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
        thisIV = IVobj(x,path,flow = 0.17, B = 10, P = 2.7e-4, V = 55,
                            I = 1.21, pos=4, ind=ind)
        # subplot(2,1,1)

        #     if ((thisIV.filter == 2) or (thisIV.filter == 2)) and (thisIV.sweepRate == 2):
        figure(1)
        isLabelonPlot = unicode(thisIV.label) in gca().get_legend_handles_labels()[1]
        thisIV.plotTrace(raw=False, deriv1 = False, deriv2 = False, label=not(isLabelonPlot))
        thisIV.plotTrace(raw=True, deriv1 = False, deriv2 = False, label=not(isLabelonPlot))
        
        figure(2)
        thisIV.plotEEDF()
        thisIV.findTe(2,8)
        print thisIV.ind, thisIV.filter, thisIV.sweepRate, thisIV.Vf, thisIV.Vp, thisIV.Te
        
        
    figure(1)
    manualVmeas = [-1.243]
    manualImeas = [.0778]
    manualV = [manualVmeas[i] - manualImeas[i] for i in range(len(manualVmeas))]
    manualI = [x/108. for x in manualImeas]
    plot(manualV, manualI, 'k*', markersize=20, label="Manual Points")
    legend(loc="best")
    text(-2, 1e-4, 'Penning, 10 G, R = 4 cm',fontsize=20)
    
    show()
   

main()
