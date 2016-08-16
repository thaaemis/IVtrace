from scipy import signal
import numpy as np
import os
from pylab import *
from IVclass import IVobj

def main():
    cwd = os.getcwd()
    
    path = "P:/departments/HTX/bkraus/ExpCampaign_EPionFlow_Aug2016/ProbeDemo_GodyakTest_081616/"
    os.chdir(path)
    files = os.listdir(path)
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
        plot(-5,-5,'g',
              -5, -5, 'b',
              -5, -5, 'r',
              -5, -5, 'c',lw=3)
        legend(['Med Filter, 1 ms', 'Med Filter, 2 ms', 'Med Filter, 4 ms', 'Med Filter, 8 ms',
                'Weak Filter, 1 ms', 'Weak Filter, 2 ms'], loc='best')
        if ((thisIV.filter == 2) or (thisIV.filter == 2)) and (thisIV.sweepRate == 2):
            figure(1)
            thisIV.plotTrace(raw=False, deriv1 = False, deriv2 = False)
            thisIV.plotTrace(raw=True, deriv1 = False, deriv2 = False)
            figure(2)
            thisIV.plotEEDF()
        print thisIV.ind, thisIV.B, thisIV.filter, thisIV.sweepRate, thisIV.Vf, thisIV.Vp, thisIV.ionSat
    try:

        show()
    except ValueError:
        print("Didn't show")
   
    os.chdir(cwd)

main()
