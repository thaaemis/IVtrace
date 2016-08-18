from scipy import signal
import numpy as np
import os
from pylab import *
from IVclass import IVobj

def main():
    cwd = os.getcwd()
    
    path = "P:/departments/HTX/bkraus/ExpCampaign_EPionFlow_Aug2016/LPEP_ThrusterIonFlow_081116/"
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
            thisIV = IVobj(x,path,flow = 2.7,B = 2,P = 6e-5,V = 250.6,
                            I = 0.616,pos = 10,ind = ind)
        else:
            thisIV = IVobj(x,path,flow = 3., B = 0, P = 6e-5, V = 68.7,
                            I = 0.73, pos=5, ind=ind)

        frame = 0 if thisIV.B == 0 else 2
        # if (thisIV.filter == 2) and (thisIV.sweepRate == 2):
        figure(frame+1)
        thisIV.plotTrace()
        figure(frame+2)
        thisIV.plotEEDF()
        thisIV.findTe(4,10)

        print thisIV.ind, thisIV.B, thisIV.filter, thisIV.sweepRate, thisIV.Vf, thisIV.Vp, thisIV.Te
    try:

        # subplot(2,1,1)
        # text(5,5e-3, 'No Field', fontsize=20)
        # subplot(2,1,2)
        # text(0,9e-3, 'Field Current = 2 A', fontsize=20)
        show()
    except ValueError:
        print("Didn't show")
   

main()
