# important features:
#   read from file DONE
#   convert file to lists DONE
#   check for current saturation DONE
#   Analyze raw curve for temperatures, Vp,Vf...? 
#   Convert to electron current DONE
#   Derivative, smoothing to get to EEDF DONE

# standard process for filtered Godyak datafile:
    # define IVobj(filename, path, flow... )
    # getFiltData(), loadFile(), correctIV()
    # getVp() and other analysis techniques (EEDF, etc.?)
    # plotTrace()

from scipy import signal
import numpy as np
from pylab import *
# from savitsky import savitzky_golay

class IVlist:
   def __init__(self):
      self.traces = []

   def sortParams(self):
      self.B,self.P,self.R,self.probe = set(),set(),set(),set()
      for x in self.traces:
         self.B.add(x.B)
         self.P.add(x.P)
         self.R.add(x.R)
         self.probe.add(x.probetype)

      self.paramDict = dict()
      for b in self.B:
         for p in self.P:
            for r in self.R:
               for probe in self.probe:
                  tmp = IVlist()
                  for x in self.traces:
                     if x.B == b and x.P == p and x.R == r:
                        tmp.traces.append(x)
                  self.paramDict[(b,p,r,probe)] = tmp
      return None

   def getProfiles(self,b,p,probe):
      R,Vp,Vf,n,Te = [],[],[],[],[]
      Vperr,Vferr,nerr,Teerr = [],[],[],[]

      for eachR in self.R:
         validTraces = self.paramDict[(b,p,eachR,probe)]
         tmpVp,tmpVf,tmpn,tmpTe = [],[],[],[]
         for x in validTraces:
            tmpVp.append(x.smoothVp)
            tmpVf.append(float(x.GodyakStats['Vf']))
            tmpn.append(float(x.GodyakStats['n (cm-3)']))
            tmpTe.append(float(x.GodyakStats['Te(eV)']))
         Vp.append(mean(tmpVp))
         Vperr.append(std(tmpVp))
         R.append(eachR)
         Vf.append(mean(tmpVf))
         Vferr.append(std(tmpVf))
         n.append(mean(tmpn))
         nerr.append(std(tmpn))
         Te.append(mean(tmpTe))
         Teerr.append(std(tmpTe))
      return R,Vp,Vperr,Vf,Vferr,n,nerr,Te,Teerr
   
class IVobj:
    def __init__(self,filename,path,flow,B,P,V,I,pos,ind):
       self.probetype = 'LP'
       self.filename = str(filename)
       self.path = str(path)
       self.flow = float(flow) # Xe gauge reading
       self.B = float(B) # in mV from gauge
       self.P = P # mTorr in chamber
       self.Vcathode = V # bias provided to cathode, V
       self.Itot = I # total plasma current
       self.R = pos
       self.ind = ind
       # for plotting:
       self.colorDict = {(2,1):'g',(2,2):'b',(2,4):'m',(2,8):'k',(4,1):'r',(4,2):'c'}

       self.getData()
       self.loadFile()
       self.correctIV()
       self.Ielectron()
       self.getVf()
       self.getVp()

    def getData(self,path = ''):
        if path == '':
            path = self.path
        with open(path+'/RawData'+self.filename,'r') as f:
            self.rawtext = f.read()
        with open(path+'/'+self.filename,'r') as f:
            self.filttext = f.read()

        return None

    def loadFile(self):
        rawTxt = self.rawtext
        filtTxt = self.filttext
        
        rawLines = rawTxt.split("\n")
        filtLines = filtTxt.split("\n")
        filtLines.pop(0)
        labels, params, V, J = [],[],[],[] # in all files
        derivJ, deriv2J, E, EEPF = [],[],[],[] # only in filtered files
        Vraw, JrawAvg = [], [] # only in raw file
        initialLine = rawLines.pop(0)
        
        rawHeadings = initialLine.split('\t')
        numTraces = rawHeadings.count('J')
        Jraw = np.zeros((numTraces, len(rawLines)-1))
        
        for line in filtLines:
            if line.count('\t') < 2:
                break
            if line.count('Parameter') < 0:
                continue
            vals = line.split('\t')
            labels.append(vals[0])
            params.append(vals[1])
            V.append(float(vals[2]))
            J.append(float(vals[3]))
            derivJ.append(float(vals[4]))
            deriv2J.append(float(vals[5]))
            E.append(vals[6])
            EEPF.append(vals[7])

        for i in range(len(rawLines)):
            line = rawLines[i]
            if line.count('\t') < 2:
                break
            vals = rawLines[i].split('\t')
            Vraw.append(float(vals[2]))
            JrawAvg.append(float(vals[3]))
            if numTraces > 0:
                for j in range(numTraces):
                    Jraw[j,i] = vals[4+j] # 0 : trace, 1 : point

        labels = list(filter(None, labels))
        params = params[:len(labels)] 
        GodyakParams = {}
        for i in range(0,len(labels)):
            GodyakParams[labels[i]] = params[i]
        self.GodyakStats = GodyakParams
        self.labels = list(filter(None,labels))
        self.params = params[:len(self.labels)]
        self.fileV = list(filter(None,V))
        if (np.shape(Jraw)[0] > 0):
            Jraw = np.mean(Jraw, 0)
        else:
            Jraw = JrawAvg
        self.Jraw = Jraw
        self.rawFileV = Vraw

        self.J = list(filter(None,J))
        self.derivJ = list(filter(None,derivJ))
        self.deriv2J = list(filter(None,deriv2J))
        self.fileE = [float(x) for x in list(filter(None,E))]
        self.fileEEPF = list(filter(None,EEPF))
        self.fileEEPF = self.fileEEPF[:len(self.fileE)]
        
        sweepRate = self.filename[self.filename.find('ms')-1] # particular to this expt
        self.filter = int(self.GodyakStats['filter '][0])
        self.sweepRate = int(sweepRate)

        return None

    def correctIV(self):
        VfGodyak = float(self.GodyakStats['Vf'])
        VpGodyak = float(self.GodyakStats['Vs'])
        self.Vraw = [float(x) + VfGodyak for x in self.rawFileV]
        self.V = [-float(x)+ VpGodyak for x in self.fileV]
        probeArea = float(self.GodyakStats['A(mm2)'])
        self.I = [x*probeArea for x in self.J]
        self.Iraw = [x*probeArea for x in self.Jraw]
        return None

    def getDerivative(self): 
        def deriv(x, y):
            x = np.asarray(x)
            y = np.asarray(y)
            return np.gradient(y)/np.gradient(x)
        self.Iprime = deriv(self.V, self.I)
        self.I2prime = deriv(self.V, self.Iprime)
      
    def Ielectron(self,current = 0): # normalize current 
        if current == 0:
            current = self.I
        try:
            self.ionSat
        except AttributeError:
            self.getIonSat()
        self.Ie = [x + self.ionSat for x in current]        
        return self.Ie

    def getVp(self):
        try:
            deriv = self.Iprime
        except AttributeError:
            self.getDerivative()
        self.Vp = self.V[find(self.Iprime == max(self.Iprime))]
        return self.Vp

    def getVf(self,smooth = False): # iffy if I doesn't cross 0
        if smooth:
            trace = np.asarray(self.smoothJ)
        else:
            trace = np.asarray(self.I)
        # establish that Vf is closest to zero point of current    
        if min(trace) < 0:
            self.Vf = self.V[find(abs(trace) == min(abs(trace)))]
        elif min(self.I) > 0:
            self.Vf = self.V[find(trace == min(trace))]

    def getIonSat(self):
        if min(self.I) < 0:
            # ion saturation is either min point, or most negatively biased point
            ionSat = abs(min(min(self.I),-abs(self.I[0])))
        elif min(self.I) > 0:
            ionSat = abs(min(self.I))
        else:
            ionSat = 0
        self.ionSat = ionSat
        return ionSat
   
   
    def plotTrace(self, deriv1=True, deriv2=True, raw=False):
        colorDict = self.colorDict
        # if self.B == 0:
        #     subplot(2,1,1)
        # else:
        #     subplot(2,1,2)
        V = self.Vraw if raw else self.V
        I = self.Iraw if raw else self.I
        markStr = ':' if raw else '-'    
        
        fig = plot(V, I,markStr, lw=3, 
            color=colorDict[(self.filter, self.sweepRate)]) #, self.V,self.smoothJ,'r',lw=3)

        # add in plasma potential as a dotted line     
        # plot([self.smoothVp,self.smoothVp],[min(self.rawJavg),max(self.rawJavg)],'-')
        yscale('log')
        if deriv1:
            plot(self.V,self.Iprime, '--',color=colorDict[(self.filter, self.sweepRate)],lw=1)
        if deriv2:
            plot(self.V,self.I2prime,'-.',color=colorDict[(self.filter, self.sweepRate)],lw=1)
        # self.plotlabel()
        ylabel('Current [A]',fontsize=20)
        xlabel('Bias Voltage [V]',fontsize=20)
        # text(self.Vp, max(self.I)/10, str(self.ind), color=colorDict[(filter, sweepRate)])
        # show()
        return None

    def plotlabel(self):
        title('B ='+str(self.B)+', P ='+str(self.P)+', R ='+str(self.R)+' mm')
        return None

    def plotEEDF(self):
        # if self.B == 0:
        #     subplot(2,1,1)
        # else:
        #     subplot(2,1,2)

        try:    
            EEDF = self.EEDF
        except AttributeError:
            self.getEEDF()

        try:
            plot(self.E,self.EEDF, color=self.colorDict[(self.filter, self.sweepRate)])
            yscale('log')
            ylabel('EEDF [arb]')
            xlabel('Energy [eV]')
        except ValueError:
            print("Didn't plot "+str(self.ind))

        return None

    def getEEDF(self):
        try:
            GodyakEEPF = self.fileEEPF
            EEDF = self.I2prime
        except AttributeError:
            EEDF = self.getDerivative

        self.E = [self.Vp - x for x in self.V]
        # for i in range(0,len(self.E)):
        #     if self.E[i] < 0:
        #         break
        # self.E = self.E[:i]
        self.EEDF = EEDF # [:i]
        return None   

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    # Old things down here, mostly for raw files
    def recalcJavg(self,recalc=True): # for raw data only
        mean,err = [],[]
        try:
            self.Jarray[0]
        except IndexError:
            self.Jerr = False
            print('Data with only one voltage sweep')
            self.rawJavg = self.rawJavg
            return None
        for ii in range(0,len(self.Jarray[0])):
            tmp = []
            for trace in self.Jarray:
                try:
                    tmp.append(trace[ii])
                except IndexError:
                    continue
            if recalc:
                mean.append(np.mean(tmp))
                err.append(np.std(tmp))
        if recalc:
            self.rawJavg = mean
        else:
            self.rawJavg = self.rawJavg
            self.rawJerr = err
        return None

    def smoothIV(self,pointrange = 100,derivative=False,order=2,EEDFfilt=50):
        # old method: smoothed = signal.medfilt(self.rawJavg,pointrange+1)
        Javg = np.asarray(self.rawJavg)
        smoothed = savitzky_golay(Javg,pointrange+1,2)
        self.smoothJ = smoothed
        print(len(smoothed),len(self.rawV))
        if derivative:
            # plot(self.rawV,self.rawJavg,self.rawV,smoothed)
            # show()
            # self.derivJ = self.difftrace(self.rawV,self.smoothJ)
            self.derivJ = savitzky_golay(self.filtJ,EEDFfilt+1,3,1)
            if order >= 2:
                # self.deriv2J = self.difftrace(self.rawV,self.derivJ)
                self.deriv2J = savitzky_golay(self.rawJavg,EEDFfilt+1,3,2)
        if derivative and order == 2:
            res = self.deriv2J
        elif derivative and order == 1:
            res = self.derivJ
        else:
            res = self.smoothJ
        return res

    def fileOK(self): # ???
        IrangetoImax = {100:100e-3,1000:10e-3,10000:1e-3} # ?? questionable
        maxcurrent = IrangetoImax[float(self.GodyakStats['I range'])]
        meascurrent = max(self.rawJavg)
        if meascurrent > 1.05*maxcurrent:
            return False
        else:
            return True
