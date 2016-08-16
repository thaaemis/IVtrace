from scipy import signal
import numpy as np
import os
from pylab import *
from IVclass import IVobj

def main():

    cwd = os.getcwd()
   
    path = "P:/departments/HTX/bkraus/ExpCampaign_EPionFlow_Aug2016/LPEP_ThrusterIonFlow_081116/"
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
        if '250.6V' in x:
            thisIV = IVobj(x,'LP',flow = 2.7,B = 2,P = 6e-5,V = 250.6,
                            I = 0.616,pos = 10,ind = ind)
        else:
            thisIV = IVobj(x,'LP',flow = 3., B = 0, P = 6e-5, V = 68.7,
                            I = 0.73, pos=5, ind=ind)
        subplot(2,1,1)
        plot(-5,-5,'g',
              -5, -5, 'b',
              -5, -5, 'm',
              -5, -5, 'k',
              -5, -5, 'r',
              -5, -5, 'c',lw=3)
        legend(['Med Filter, 1 ms', 'Med Filter, 2 ms', 'Med Filter, 4 ms', 'Med Filter, 8 ms',
                'Weak Filter, 1 ms', 'Weak Filter, 2 ms'], loc='best')

        
        thisIV.getfiltdata(path)
        # thisIV.getrawdata(path)
        # thisIV.readfile(True)
        thisIV.readfile(False)
        thisIV.correctVaxis()
        thisIV.getVp()
        thisIV.plotTrace(False, False)
        subplot(2,1,1)
        text(0,5e-3, 'No Field', fontsize=20)
        subplot(2,1,2)
        text(0,9e-3, 'Field Current = 2 A', fontsize=20)

        
        print thisIV.B, thisIV.filter, thisIV.sweepRate, thisIV.Vp
        
    show()
   
   
   
   
    os.chdir(cwd)

main()
