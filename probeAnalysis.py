import re
import numpy as np
import os, sys
from scipy import signal

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
class Probes(object):
    def __init__(self):
        self.time = None
        self.p = None
        self.T1 = None
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
        
    def getProbe(self,n):
        if n in self.list:
            index=self.list.index(n)
        else:
            index = 0
            print("probe is not in the list, reporting probe #"+ str(self.list[index]))
        #print("Index========================",index)
        
        vars= Probes()
        
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
class processFileReturn(object):
    def __init__(self,data,probeList,probeLoc,probeLoctn,times, varName):
        self.data = data
        self.probeList = probeList
        self.probeLoc = probeLoc
        self.probeLoctn = probeLoctn
        self.times = times
        self.varName = varName

    
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
    return(processFileReturn(data,probeList,probeLoc,probeLoctn,times,varName))
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
    ddt=[]
    ddt.append((signal[1]-signal[0])/(time[1]-time[0]))     #handling first data point derivate (first order-forward)
    for i in range(1,len(signal)-1):
        ddt.append((signal[i+1]-signal[i-1])/(time[i+1]-time[i-1]))
    ddt.append((signal[-1]-signal[-2])/(time[-1]-time[-2])) #handling last data point derivate (first order-backward)
    return(np.asarray(ddt))
    
    
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
  filtered_audio = signal.filtfilt(filter_coefs, [1], y)
  return filtered_audio
  
  
### select probes with specific criteria (e.q. locaiton)
def selectProbes(loctn, s, walldist):
  index=[]
  for j in range(loctn.shape[0]):
        if (loctn[j][0]>s[0] and loctn[j][0]<s[1] and loctn[j][1]>walldist[0] and loctn[j][1]<walldist[1]):
            index.append(j)
            
  print(index)
  print(loctn[index])
  return index
  
class Gamma(object):
    def __init__(self,s,wallDist,gamma):
        self.s = s
        self.wallDist = wallDist
        self.gamma = gamma
        

### cluster probes corresponding to their locations, outputs an array with corresponding bin numbers, last col remains gamma
def clusterProbes(gamma, thresh):
    print(gamma.shape)
    gammaIndex=np.insert(gamma,[0,0],np.ones([len(gamma[:,0]),2]),axis=1)
    print(gammaIndex.shape)
    
    ####inserted two columns to the gamma at the beginning, to preserve s and w to be used later, first column is the s _index, the second column would be the w index 
    # thresh 0 - s, thresh 1- w, thresh 3 - z
    #sorting based on s
    gammaIndex=gammaIndex[np.argsort(gammaIndex[:,2])]
    
    print('gammaIndex -sorted s\n',gammaIndex.round(6))
    #calculate the diff in sorted s1
    diffs=np.diff(gammaIndex[:,2])
    diffIndex=np.argwhere(diffs>thresh[0]).flatten()
    diffIndex = np.hstack([0, diffIndex + 1, len(gammaIndex)+1])
    counters=0
    for a, b in zip(diffIndex[:-1], diffIndex[1:]):  # iterate over index pairs
        print('ab=',a,',',b)
        gammaIndex[a:b,0]=counters
        counters=counters+1
        
        subgammaIndex=gammaIndex[a:b,:]
        subgammaIndex=subgammaIndex[np.argsort(subgammaIndex[:,3])]
        gammaIndex[a:b,:]=subgammaIndex
        diffw=np.diff(subgammaIndex[:,3])
        diffwIndex=np.argwhere(diffw>thresh[1]).flatten()
        diffwIndex = np.hstack([0, diffwIndex + 1, len(subgammaIndex)])
        print('sub gamma\n',subgammaIndex.round(6))
        counterw=0
        for c, d in zip(diffwIndex[:-1], diffwIndex[1:]):  # iterate over index pairs
            print('cd=',c,',',d)
            gammaIndex[c+a:d+a,1]=counterw
            counterw=counterw+1
        print('updated gammaIndex in sub loop\n',gammaIndex.round(6))
    GammaList=[]
    print(np.max(gammaIndex[:,0]))
    for i in range(int(np.max(gammaIndex[:,0]))+1):
        #filter this location
        gammas=gammaIndex[gammaIndex[:,0]==i,:]
        #print('gammas\n',gammas.round(3))
        s=np.mean(gammas[:,2])
        
        #print(s)
        w=[] # list of wall distances  for this s
        g=[] # list of gamma for this s
        for j in range(int(np.max(gammas[:,1]))+1):
            gammaw=gammas[gammas[:,1]==j,:]
            #print(gammaw.round(4))
            w.append(np.mean(gammaw[:,3]))
            g.append(np.mean(gammaw[:,5]))
        GammaList.append(Gamma(s,w,g))
        # print(GammaList[i].s)
        # print(GammaList[i].wallDist)
        # print(GammaList[i].gamma)
    
    return GammaList
