import re
import numpy as np
import os, sys
from scipy import signal
from skimage import filters
import matplotlib.pyplot as pl
from scipy.ndimage.filters import generic_filter

def yes_or_no(question):

    reply = str(input(question+' (y/n): ')).lower().strip()

    if reply[0] == 'y':

        return True

    if reply[0] == 'n':

        return False

    else:

        return yes_or_no("Uhhhh... please enter ")
## transfers the coordinates
## calculates s,wallDist and z-axiz
## calculates Rotation tensor which can be used if needed! now most of the parameters are already rotated. 
class coordinateTransfer(object):
    def __init__(self):
        #print(os.path.dirname(os.path.realpath(__file__)))
        self.xc, self.yc, self.s, self.slope=np.loadtxt(os.path.dirname(os.path.realpath(__file__)) + "/curve_py",unpack=True)
    def transferAll(self, locs):
        probeLoctn=np.empty((0,3), int)
        for j in range(locs.shape[0]):
            probeLoctn=np.vstack((probeLoctn,self.transfer(locs[j])))
        return(probeLoctn)
    def transfer(self,loc):
        x=loc[0]
        y=loc[1]
        z=loc[2]
        distMin1=1e10
        distMin2=1e10
        i1=0
        i2=0
        unitx=np.array([1, 0, 0])
        Rx=np.array(([1,0,0],[0,-1,0],[0,0,-1]))
        for i in range(len(self.xc)):
                dist = np.sqrt(pow((self.xc[i]-x),2)+pow((self.yc[i]-y),2))
                if (dist < distMin1):
                        distMin2=distMin1
                        i2 = i1        
                        distMin1=dist
                        i1 = i
        if (i1==0):
                i2 = 1
        ## first point on the curve
        x1 = self.xc[i1]
        y1 = self.yc[i1]
        ##second point on the curve
        x2 = self.xc[i2]
        y2 = self.yc[i2]
        ## slopes associated with those two points
        m1 = self.slope[i1]
        m2 = self.slope[i2]
        ## surface distance associated with two points
        s1 = self.s[i1]
        s2 = self.s[i2]
        ## interpolating the point and slope
        m = (y2-y1)/(x2-x1)
        b = 1
        a = -m
        c = -y1+m*x1
        xi = (b*( b*x -a*y)-a*c)/(pow(a,2)+pow(b,2))
        yi = (a*(-b*x +a*y)-b*c)/(pow(a,2)+pow(b,2))
        wallDist = np.absolute(a*x+b*y+c)/np.sqrt((pow(a,2)+pow(b,2)))
        mi = m1+ (xi-x1)*(m2-m1)/(x2-x1)
        si = s1+ (xi-x1)*(s2-s1)/(x2-x1)
        tv = np.array([1,mi,0])
        tvnorm = tv / np.linalg.norm(tv)
        if (mi < 0):
                sign=-1
        else:
                sign=1
        Rz= np.array(([np.dot(unitx,tvnorm),sign*np.linalg.norm(np.cross(tvnorm,unitx)),0],[-sign*np.linalg.norm(np.cross(tvnorm,unitx)),np.dot(unitx,tvnorm),0],[0,0,1]));
        
        if (pow(xi+0.22,2)+pow(yi+0.4489988,2) < 0.25): #inside the circle--> pressure side, see the coordinate excel file
                self.R=np.dot(Rx,Rz);        #first rotate around z-axiz, then rotate around x-axis, order is important!
        else: #on suction surface, no need to rotate around x
                self.R=Rz;

        return(np.array([si, wallDist,z]))

# creates an object that contains all probe data
class ProbeVars(object):
    def __init__(self):
        self.time = None
        self.p = None
        self.T1 = None
        self.T1mean = None
        self.T1prime = None
        self.T2 = None
        self.Utn = None
        self.U = None
        self.D = None
        self.filterDA = None
        self.filterDG = None
        self.filterQA = None
        self.gradU = None
        self.Q = None
        self.vorticity = None
        self.list = None
        self.loc = None
        self.loctn = None
        self.u = None
        self.v = None
        self.w = None
        self.uprime = None
        self.vprime = None
        self.wprime = None
        self.umean = None
        self.vmean = None
        self.wmean = None
        self.detector = None
        
    def getProbe(self,n):
        ## to get probe number, deprecated
        # if n in self.list:
        #     index=self.list.index(n)
        # else:
        #     index = 0
        #     print("probe is not in the list, reporting probe #"+ str(self.list[index]))
        #print("Index========================",index)
        index = n
        vars= ProbeVars()
        
        vars.p = self.p[index,:]
        vars.T1 = self.T1[index,:]
        vars.T2 = self.T2[index,:]
        vars.U = self.U[index,:,:]
        vars.Utn = self.Utn[index,:,:]
        vars.D = self.D[index,:]
        vars.filterDA = self.filterDA[index,:]
        vars.filterDG = self.filterDG[index,:]
        vars.filterQA = self.filterQA[index,:]
        vars.gradU = self.gradU[index,:,:]
        vars.Q = self.filterQA[index,:]
        vars.vorticity = self.vorticity[index,:]
        
        vars.list = self.list[index]
        vars.loc = self.loc[index]
        vars.loctn = self.loctn[index]
        vars.time=self.time
        
        return(vars)
        

#######################################
############ Processing files #########
class Probe(object):
    def __init__(self):
        self.data = None
        self.list = None
        self.loc = None
        self.loctn = None
        self.time = None
        self.varName = None
    def getProbe(self,n):
            ## to get probe number, deprecated
            # if n in self.list:
            #     index=self.list.index(n)
            # else:
            #     index = 0
            #     print("probe is not in the list, reporting probe #"+ str(self.list[index]))
            #print("Index========================",index)
            index=n
            probe= Probe()
            probe.list = self.list[index]
            probe.loc = self.loc[index]
            probe.loctn = self.loctn[index]
            probe.time=self.time
            if self.data.ndim == 2:
                probe.data = self.data[index,:]
            else:
                probe.data = self.data[index,:,:]
            
            return(probe)

    def append(self,probe):
        if self.loc is None:
            self.data =probe.data
            self.list = probe.list
            self.loc = probe.loc
            self.loctn =probe.loctn
            self.time = probe.time
            self.varName = probe.varName
        else:
            self.data = np.concatenate((self.data,probe.data))
            self.list = np.concatenate((self.list,probe.list))
            self.loc = np.concatenate((self.loc,probe.loc))
            self.loctn = np.concatenate((self.loctn,probe.loctn))


    
def processFile(probeInfo,probeName,initTime,varName,nComp):
    file = "postProcessing/" + probeName + "/" + initTime + "/" + varName
    print("started reading " + file)
    probeList=[]
    dataTemp=[]
    xp=[]
    yp=[]
    zp=[]
    #probeList.append(probeInfo[0])  # appending first probe ORG
    probeList=[]
    colList=[]
    AvailableProbeCounter=-1
    probeListInput=range(probeInfo[0],probeInfo[1]+1,probeInfo[2])
    probesFound=0
    #with open(file) as f:
    #    head = [next(f) for x in range(5000)]
    #for line in head:
    with open(file) as f:
        for line in f:
            if(line[0]!='#'):
                row = [y for y in re.split('[(|)]|\s+',line) if y]
                dataTemp.append(row[:])
            else: #if (not probesFound):   #finding the probe location based on probe number
                pointdata=[y for y in re.split('[(|)]|\s+',line) if y]
                #if ((len(pointdata)> 4) and  (float(pointdata[2])==probeList[-1])): ##ORG
                if (len(pointdata)> 4 and len(pointdata)< 7 ):
                    #print(pointdata)
                    AvailableProbeCounter=AvailableProbeCounter+1
                    if (float(pointdata[2]) in probeListInput):
                        xp.append(float(pointdata[3]))
                        yp.append(float(pointdata[4]))
                        zp.append(float(pointdata[5]))
                        probeList.append(float(pointdata[2]))
                        colList.append(AvailableProbeCounter)
                    # if (probeList[-1]<probeInfo[1]):
                            # probeList.append(float(pointdata[2])+probeInfo[2])
                    # else:
                            # probesFound = 1
    #xp0=[(x-xp[0]) for x in xp]

    #print(probeList)
    probeLoc=np.array([list(x) for x in zip(xp, yp, zp)])
    coordTransfer = coordinateTransfer()    #initializing coordinate transfer
    probeLoctn = coordTransfer.transferAll(probeLoc)
    dataZipped=list(zip(*dataTemp))
    times = np.array(dataZipped[0],dtype=float)
    #data=np.array(list(dataZipped[int(x)*nComp+1] for x in probeList),dtype=float)
    data=np.array(list(dataZipped[int(x)*nComp+1] for x in colList),dtype=float)
    for comp in range(1,nComp):
       data=np.dstack((data, np.array(list(dataZipped[int(x)*nComp+comp+1] for x in colList),dtype=float)))
    #data = np.dstack((np.array(list(dataZipped[int(x)*3+1] for x in probeList),dtype=float), 
    #       np.array(list(dataZipped[int(x)*3+2] for x in probeList),dtype=float), 
    #       np.array(list(dataZipped[int(x)*3+3] for x in probeList),dtype=float)))
    print("finished reading " + file)
    probe=Probe()
    probe.data=data
    probe.time=times
    probe.list=probeList
    probe.loc=probeLoc
    probe.loctn=probeLoctn
    probe.varName=varName
    return(probe)
#
### Data structure 
###           times     
###         -------> 
###        |
### probes |
###        ˅
### i.e,
### probes are in the rows (first index data[#,:,:])
### times are in the columns(second index data[:,#,:])
### components are in the third dimension (data[:,:,#])
### e.g. data[:,:,1] for velocity gives the U-component for all saved probes in time
###########

def ddt(signal,time):
    ddt = np.gradient(signal, time,axis=1)
    # if signal.ndim == 1:
    #     ddt=[]
    #     ddt.append((signal[1]-signal[0])/(time[1]-time[0]))     #handling first data point derivate (first order-forward)
    #     for i in range(1,len(signal)-1):
    #         ddt.append((signal[i+1]-signal[i-1])/(time[i+1]-time[i-1]))
    #     ddt.append((signal[-1]-signal[-2])/(time[-1]-time[-2])) #handling last data point derivate (first order-backward)
    # elif signal.ndim == 2:
        

        # ddt=np.zeros(signal.shape)
        # ddt[:,1]=(signal[:,1]-signal[:,0])/(time[1]-time[0])     #handling first data point derivate (first order-forward)
        # for i in range(1,signal.shape[1]-1):
        #     ddt=np.append(ddt,(signal[:,i+1]-signal[:,i-1])/(time[i+1]-time[i-1]),axis=1)
        # ddt=np.append(ddt,(signal[:,-1]-signal[:,-2])/(time[-1]-time[-2]),axis=1) #handling last data point derivate (first order-backward)
    
    return(ddt)
    

    
#################
##### high pass filter
def highpass_filter(y, sr,filter_stop_freq,filter_pass_freq,filter_order):
  # filter_stop_freq = 2300  # Hz
  # filter_pass_freq = 2510  # Hz
  # filter_order = 1001

  # High-pass filter
  nyquist_rate = sr / 2.
  desired = (0, 0, 1, 1)
  bands = (0, filter_stop_freq, filter_pass_freq, nyquist_rate)
  filter_coefs = signal.firls(filter_order, bands, desired, nyq=nyquist_rate)

  # Apply high-pass filter
  filtered_audio = signal.filtfilt(filter_coefs, [1], y,axis=0)
  
  return filtered_audio
  
  
### select probes with specific criteria (e.q. locaiton)
def selectProbes(loctn, s, walldist):
  index=[]
  for j in range(loctn.shape[0]):
        if (loctn[j][0]>s[0] and loctn[j][0]<s[1] and loctn[j][1]>walldist[0] and loctn[j][1]<walldist[1]):
            index.append(j)
            
  #print(index)
  #print(loctn[index])
  return index
  ### select probes with specific criteria (e.q. locaiton)
def selectProbeIndexes(probeIndex,s,walldist):
  index=[]
  for j in range(probeIndex.shape[0]):
        if (probeIndex[j][0]>=s[0] and probeIndex[j][0]<=s[1] and probeIndex[j][1]>=walldist[0] and probeIndex[j][1]<=walldist[1]):
            index.append(j)
            
  #print(index)
  #print(loctn[index])
  return index
class Gamma(object):
    def __init__(self,s,wallDist,gamma,delta99):
        self.s = s
        self.wallDist = wallDist
        self.gamma = gamma
        self.delta99 = delta99
        

### cluster probes corresponding to their locations, outputs an array with corresponding bin numbers, last col remains gamma
def clusterProbes(gamma, thresh):
    #print(gamma.shape)
    sref, delta99ref=np.loadtxt(os.path.dirname(os.path.realpath(__file__)) + "/delta99",unpack=True)
    gammaIndex=np.insert(gamma,[0,0],np.ones([len(gamma[:,0]),2]),axis=1)
    #print(gammaIndex.shape)
    
    ####inserted two columns to the gamma at the beginning, to preserve s and w to be used later, first column is the s _index, the second column would be the w index 
    # thresh 0 - s, thresh 1- w, thresh 3 - z
    #sorting based on s
    gammaIndex=gammaIndex[np.argsort(gammaIndex[:,2])]
    
    #print('gammaIndex -sorted s\n',gammaIndex.round(6))
    #calculate the diff in sorted s1
    diffs=np.diff(gammaIndex[:,2])
    diffIndex=np.argwhere(diffs>thresh[0]).flatten()
    diffIndex = np.hstack([0, diffIndex + 1, len(gammaIndex)+1])
    counters=0
    for a, b in zip(diffIndex[:-1], diffIndex[1:]):  # iterate over index pairs
        #print('ab=',a,',',b)
        gammaIndex[a:b,0]=counters
        counters=counters+1
        
        subgammaIndex=gammaIndex[a:b,:]
        subgammaIndex=subgammaIndex[np.argsort(subgammaIndex[:,3])]
        gammaIndex[a:b,:]=subgammaIndex
        diffw=np.diff(subgammaIndex[:,3])
        diffwIndex=np.argwhere(diffw>thresh[1]).flatten()
        diffwIndex = np.hstack([0, diffwIndex + 1, len(subgammaIndex)])
        #print('sub gamma\n',subgammaIndex.round(6))
        counterw=0
        for c, d in zip(diffwIndex[:-1], diffwIndex[1:]):  # iterate over index pairs
            #print('cd=',c,',',d)
            gammaIndex[c+a:d+a,1]=counterw
            counterw=counterw+1
        #print('updated gammaIndex in sub loop\n',gammaIndex.round(6))
    GammaList=[]
    gammaarray=np.empty( shape=(0, 4) )

    #print(np.max(gammaIndex[:,0]))
    for i in range(int(np.max(gammaIndex[:,0]))+1):
        #filter this location
        gammas=gammaIndex[gammaIndex[:,0]==i,:]
        #print('gammas\n',gammas.round(3))
        s=np.mean(gammas[:,2])
        #delta99=0.0007299381*np.exp(4.025940745*s)
        delta99 = np.interp(s,sref,delta99ref)
        #print(s)
        w=[] # list of wall distances  for this s
        g=[] # list of gamma for this s
        for j in range(int(np.max(gammas[:,1]))+1):
            gammaw=gammas[gammas[:,1]==j,:]
            #print(gammaw.round(4))
            w.append(np.mean(gammaw[:,3]))
            g.append(np.mean(gammaw[:,5]))
        GammaList.append(Gamma(s,w,g,delta99))
        warray=np.asarray(w)
        garray=np.asarray(g)
        sarray= np.ones((len(w),1))*s
        delta99array=np.ones((len(w),1))*delta99
        temp=np.column_stack( (sarray,warray,garray,delta99array))
        gammaarray=np.vstack((gammaarray,temp))
        # print(GammaList[i].s)
        # print(GammaList[i].wallDist)
        # print(GammaList[i].gamma)
    print(w)
    return GammaList, np.asarray(gammaarray)
def ProbeIndex(gamma, thresh):
    #print(gamma.shape)
    gammaIndex=np.insert(gamma,[0,0],np.ones([len(gamma[:,0]),2]),axis=1)
    #print(gammaIndex.shape)
    
    ####inserted two columns to the gamma at the beginning, to preserve s and w to be used later, first column is the s _index, the second column would be the w index 
    # thresh 0 - s, thresh 1- w, thresh 3 - z
    #sorting based on s
    gammaIndex=gammaIndex[np.argsort(gammaIndex[:,2])]
    
    #print('gammaIndex -sorted s\n',gammaIndex.round(6))
    #calculate the diff in sorted s1
    diffs=np.diff(gammaIndex[:,2])
    diffIndex=np.argwhere(diffs>thresh[0]).flatten()
    diffIndex = np.hstack([0, diffIndex + 1, len(gammaIndex)+1])
    counters=0
    for a, b in zip(diffIndex[:-1], diffIndex[1:]):  # iterate over index pairs
        #print('ab=',a,',',b)
        gammaIndex[a:b,0]=counters
        counters=counters+1
        
        subgammaIndex=gammaIndex[a:b,:]
        subgammaIndex=subgammaIndex[np.argsort(subgammaIndex[:,3])]
        gammaIndex[a:b,:]=subgammaIndex
        diffw=np.diff(subgammaIndex[:,3])
        diffwIndex=np.argwhere(diffw>thresh[1]).flatten()
        diffwIndex = np.hstack([0, diffwIndex + 1, len(subgammaIndex)])
        #print('sub gamma\n',subgammaIndex.round(6))
        counterw=0
        for c, d in zip(diffwIndex[:-1], diffwIndex[1:]):  # iterate over index pairs
            #print('cd=',c,',',d)
            gammaIndex[c+a:d+a,1]=counterw
            counterw=counterw+1
        #print('updated gammaIndex in sub loop\n',gammaIndex.round(6))
    return gammaIndex

def lt(detector,method,input,smoothingLength):
    N=smoothingLength #100
    dSmooth=signal.convolve(detector, np.ones((N,))/N, mode='valid')
    #dSmooth=signal.convolve(dSmooth, np.ones((N,))/N, mode='valid')
    #dSmooth2=running_mean(detector,N)
    dt=1e-5
    #dSmooth2=highpass_filter(detector,1/dt,1/(dt*50),2/(dt*50),1001)
    # dSmooth = generic_filter(detector, np.std, size=N)

    ##std filter
    # blur=signal.convolve(detector, np.ones((N,))/N, mode='valid')
    # blur2=signal.convolve(np.square(detector), np.ones((N,))/N, mode='valid')
    # dSmooth=np.sqrt(blur2-np.square(blur))

    # pl.figure(10)
    # pl.plot(detector,'k')
    # pl.plot(dSmooth,'-r')
    # pl.plot(dSmooth2,'-b')
    # pl.show()
    ## Laminar/Turbulent Discrimination
    if method=='otsuloc':
        ### Otsu
        thresh_otsu= filters.threshold_otsu(dSmooth)
        lt_otsu=np.where(dSmooth>thresh_otsu, 1, 0) # Laminar_Turbulent, Ostu
        lamtu = lt_otsu
        print("Otsu Threshhold={:.2e}".format(thresh_otsu))
        #pl.plot(time[int(N/2):int(-N/2+1)],lt_otsu,'-r',label="LT")
    elif method=='adaptive':
        # ### Adaptive
        block_size = 5001
        print(np.tile(dSmooth,(2,1)).shape)
        lt_adaptive = filters.threshold_adaptive(np.tile(dSmooth,(2,1)), block_size, offset=100)*1
        lamtu = lt_adaptive[0,:]
        print(lt_adaptive[0,:].shape)
        #pl.plot(time[int(N/2):int(-N/2+1)],lt,'-k',label="LT") 
    elif method=='percentageMax':
        ### Arbitrary fraction of Max D
        C = input
        
        lt_cmax=np.where(dSmooth>C*np.max(dSmooth),1,0)#0.01*24473540.92065062,1,0)#>0.01*269112747.3345178,1,0) ##6419964748.39, 1, 0) # Laminar_Turbulent, C*max(D)
        lamtu = lt_cmax
        #pl.plot(time[int(N/2):int(-N/2+1)],lt_cmax,'-b',label="LT")
    elif method=='threshglob':
        thresh = input 
        ### threshhold, can be c*Max(D) or otsu global   
        #lt_cmax=np.where(DSmooth>C*np.max(DSmooth), 1, 0) # Laminar_Turbulent, C*max(D)
        lt_cmax=np.where(dSmooth>thresh,1,0)#0.01*24473540.92065062,1,0)#>0.01*269112747.3345178,1,0) ##6419964748.39, 1, 0) # Laminar_Turbulent, C*max(D)
        lamtu = lt_cmax
        #pl.plot(time[int(N/2):int(-N/2+1)],lt_cmax,'-b',label="LT")
    return lamtu

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def general_gamma(x,xs,xe):
    return 1-np.exp(-5* np.power( (x-xs) /(xe-xs),3  ))