# important features:
#   read from file DONE
#   convert file to lists DONE
#   check for current saturation DONE
#   Analyze raw curve for temperatures, Vp,Vf...? 
#   Convert to electron current DONE
#   Derivative, smoothing to get to EEDF DONE

from scipy import signal
import numpy
from pylab import *
from savitsky import savitzky_golay

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
   def __init__(self,filename,probetype,flow,B,P,V,I,pos,ind):
      probeDict = {'EP':0,'LP':1,'RP':2,'IP':3,'BP':4}
      def LPr(r):
         return (4-r)*10 # mm to cm
      def RPr(r):
         return (4-r-4.5)*10 # mm to cm
      def reg(r):
         return r
      radfunc = [reg,LPr,RPr,reg,reg]
      self.probetype = probeDict[probetype]
      proberad = radfunc[self.probetype](pos)
      self.filename = str(filename)
      self.flow = float(flow) # Xe gauge reading
      self.B = float(B) # in mV from gauge
      self.P = P # mTorr in chamber
      self.Vcathode = V # bias provided to cathode, V
      self.Itot = I # total plasma current
      self.R = proberad
      self.ind = ind

   def getrawdata(self,path):
      with open(path+'/RawData'+self.filename,'r') as f:
         self.rawtext = f.read()
      return None

   def getfiltdata(self,path):
      with open(path+'/'+self.filename,'r') as f:
         self.filttext = f.read()
      return None

   def correctVaxis(self):
      VpGodyak = float(self.GodyakStats['Vs'])
      self.V = [-float(x)+VpGodyak for x in self.filtV]
      
      probeArea = float(self.GodyakStats['A(mm2)'])
      self.I = [x*probeArea for x in self.filtJ]
      
      return None

   def readfile(self,raw=True):
      try:    
         txt = self.rawtext if raw else self.filttext
      except AttributeError:
         print('First get data file with self.getrawdata(), self.getfiltdata()!')

      lines = txt.split("\n")
      if raw:
         labels,params,V,Javg,J = [],[],[],[],[]
         for i in range(1,len(lines)):
            line = lines[i]
            if line.count('\t') < 2:
               break
            if line.count('J') > 5:
               continue
            vals = line.split('\t')
            labels.append(vals[0])
            params.append(vals[1])
            print(vals)
            V.append(float(vals[2]))
            Javg.append(float(vals[3]))
            for i in range(4,len(vals)):
               J.append([])
               J[i-4].append(float(vals[i]))

         labels = list(filter(None,labels))
         params = params[:len(labels)] 
         GodyakParams = {}
         for i in range(0,len(labels)):
            GodyakParams[labels[i]] = params[i]
         self.GodyakStats = GodyakParams
         self.rawV = V
         self.rawJavg = Javg
         self.Jarray = J
         print(len(self.rawJavg))
         return None
      else:
         self.filtlabels,self.filtparams,self.filtV,self.filtJ = [],[],[],[]
         self.filtderivJ,self.filtderiv2J,self.filtE,self.filtEEPF = [],[],[],[]
         for i in range(1,len(lines)):
            line = lines[i]
            if line.count('\t') < 2:
               break
            if line.count('Parameter') < 0:
               continue
            vals = line.split('\t')
            self.filtlabels.append(vals[0])
            self.filtparams.append(vals[1])
            self.filtV.append(float(vals[2]))
            self.filtJ.append(float(vals[3]))
            self.filtderivJ.append(float(vals[4]))
            self.filtderiv2J.append(float(vals[5]))
            self.filtE.append(vals[6])
            self.filtEEPF.append(vals[7])
         labels = list(filter(None,self.filtlabels))
         params = self.filtparams[:len(labels)] 
         GodyakParams = {}
         for i in range(0,len(labels)):
            GodyakParams[labels[i]] = params[i]
         self.GodyakStats = GodyakParams

         self.filtlabels = list(filter(None,self.filtlabels))
         self.filtparams = self.filtparams[:len(self.filtlabels)]
         self.filtV = list(filter(None,self.filtV))
         self.filtJ = list(filter(None,self.filtJ))
         self.filtderivJ = list(filter(None,self.filtderivJ))
         self.filtderiv2J = list(filter(None,self.filtderiv2J))
         self.filtE = list(filter(None,self.filtE))
         self.filtE = [float(x) for x in self.filtE]
         self.filtEEPF = list(filter(None,self.filtEEPF))
         self.filtEEPF = self.filtEEPF[:len(self.filtE)]
         return None

   def fileOK(self):
      IrangetoImax = {100:100e-3,1000:10e-3,10000:1e-3} # ?? questionable
      maxcurrent = IrangetoImax[float(self.GodyakStats['I range'])]
      meascurrent = max(self.rawJavg)
      if meascurrent > 1.05*maxcurrent:
         return False
      else:
         return True

   def difftrace(self,x,y):
      deriv = []
      dx = x[3]-x[2]
      for i in range(1,len(y)-1):
         deriv.append((y[i+1]-y[i-1])/(2*dx))
      deriv[0] = 2*deriv[1]-deriv[2]
      deriv[-1] = 2*deriv[-2] - deriv[-3]
      return deriv

   def getDerivative(self): 
       def deriv(x, y):
           x = np.asarray(x)
           y = np.asarray(y)
           return np.gradient(y)/np.gradient(x)
       self.Iprime = deriv(self.V, self.I)
      
   def recalcJavg(self,recalc=True):
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
            mean.append(numpy.mean(tmp))
            err.append(numpy.std(tmp))
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
      
      
      
   def Ielectron(self,current):
      try:
         self.Ie = current + self.ionSat
      except AttributeError:
         self.getIonSat()
         self.Ie = current + self.ionSat
      return self.Ie

   def getVp(self):
      try:
         deriv = self.Iprime
      except AttributeError:
         # self.smoothIV(100,True,1)
         self.getDerivative()
      self.Vp = self.V[find(self.Iprime == max(self.Iprime))]
      return self.Vp

   def Vfraw(self,smooth = False):
      if smooth:
         trace = self.smoothJ
      else:
         trace = self.rawJavg
      # establish that Vf is closest to zero point of current    
      if min(trace) < 0:
         self.Vf = self.rawV[find(abs(trace) == min(abs(trace)))]
      elif min(self.rawJavg) > 0:
         self.Vf = self.rawV[find(trace == min(trace))]

   def getIonSat(self):
      if min(self.rawJavg) < 0:
         # ion saturation is either min point, or most negatively biased point
         ionSat = abs(min(min(self.rawJavg),-abs(self.rawJavg[0])))
      elif min(self.rawJavg) > 0:
         ionSat = abs(min(self.rawJavg))
      else:
         ionSat = 0
         self.ionSat = ionSat
      return ionSat
   
   def plotTrace(self,deriv1=True,deriv2=True):
      
      sweepRate = self.filename[self.filename.find('ms')-1]
      filter, sweepRate = int(self.GodyakStats['filter '][0]), int(sweepRate)
      self.filter = filter
      self.sweepRate = sweepRate

      colorDict = {(2,1):'g',(2,2):'b',(2,4):'m',(2,8):'k',(4,1):'r',(4,2):'c'}
      
   
      if self.B == 0:
         subplot(2,1,1)
      else:
         subplot(2,1,2)
         
      fig = plot(self.V,self.I,'-', lw=3, color=colorDict[(filter, sweepRate)]) #, self.V,self.smoothJ,'r',lw=3)

      # add in plasma potential as a dotted line     
      # plot([self.smoothVp,self.smoothVp],[min(self.rawJavg),max(self.rawJavg)],'-')
      yscale('log')
      if deriv1:
         ratio1 = max(self.rawJavg)/max(self.derivJ)/8
         plot(self.V,self.derivJ*ratio1,'m',lw=3)
      if deriv2:
         ratio2 = max(self.rawJavg)/max(self.deriv2J)/50
         plot(self.V,self.deriv2J*ratio2,'g',lw=4)
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
      try:    
         plot(self.E,self.EEDF)
      except AttributeError:
         self.getEEDF()
      yscale('log')
      self.plotlabel()
      show()
      return None

   def getEEDF(self):
      try:
         EEDF = self.deriv2J
      except AttributeError:
         EEDF = self.smoothIV(100,True,2,50)

      self.E = [self.smoothVp - x for x in self.V]
      for i in range(0,len(self.E)):
         if self.E[i] < 0:
            break
      self.E = self.E[:i]    
      self.EEDF = EEDF[:i]
      return None   

