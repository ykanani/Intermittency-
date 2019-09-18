# reads the data from multiple folders
import os, sys
#import statistics
import math
import matplotlib.pyplot as pl
import numpy as np
import csv
import pickle
#import cv2 as cv
from skimage import filters
import re
from scipy import signal
from scipy import stats
from scipy.interpolate import spline
from decimal import Decimal
import datetime
#os.chdir("D:\Dropbox\PhD\Python")
#import easygui
from matplotlib import rc
from probeAnalysis import ProbeVars
from probeAnalysis import Gamma
from probeAnalysis import selectProbes
from probeAnalysis import selectProbeIndexes
from probeAnalysis import processFile
from probeAnalysis import clusterProbes
from probeAnalysis import highpass_filter
from probeAnalysis import ProbeIndex
from probeAnalysis import Probe
from probeAnalysis import lt
from probeAnalysis import ddt
from probeAnalysis import yes_or_no
from probeAnalysis import general_gamma
from probeAnalysis import filterSignal
from scipy.interpolate import interp1d
from scipy.ndimage.filters import generic_filter
from scipy.optimize import curve_fit
import matplotlib as mpl; print(mpl.font_manager.get_cachedir())

#####GOOD FONTS!!!! remove comment
rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':20})
rc('text', usetex=True)
os.chdir(os.path.dirname(os.path.realpath(__file__)))
print(os.getcwd())

os.chdir('I:/PostProcess/ccc2')
os.listdir()
print(os.getcwd())
date_time = datetime.datetime.now()
print("GammaList"+str(date_time)) 
##input parameters
nu = 1.55e-5
#nu =7.7873794212218649517684887459807e-6 #for scaled length and velocity
L=0.4976 #chord
Uin=4
T=L/Uin
dt=1e-5

#####################
##################### Main
#####################


ndelta=21  # number of probes at each station (wall normal)
ns=90      # number of stations
probeInfo=[0,ndelta*ns-1,1] # [first,last,interval]
zlocs=['-0.05', '-0.04', '-0.03', '-0.02', '-0.01', '0.01', '0.02', '0.03', '0.04', '0.05']
gamma = np.empty([0,4]) # culumns: s,w,z,gamma

DAll=Probe()
UtnAll=Probe()
T1All=Probe()
filterDAAll=Probe()
for k in range(len(zlocs)):
    probeName=  "p_ss_BL_z" + zlocs[k]# "test"  # 
    print('reading data file',"D" + "_" + probeName)
    with open("postProcessing/" + "D" + "_" + probeName, "rb") as fp:   #unPickling
        DAll.append(pickle.load(fp))
    print('reading data file',"filterDA" + "_" + probeName)
    with open("postProcessing/" + "filterDA" + "_" + probeName, "rb") as fp:   #unPickling
        filterDAAll.append(pickle.load(fp))
    print('reading data file',"Utn" + "_" + probeName)
    with open("postProcessing/" + "Utn" + "_" + probeName, "rb") as fp:   #unPickling
        UtnAll.append(pickle.load(fp))
    print('reading data file',"T1" + "_" + probeName)
    with open("postProcessing/" + "T1" + "_" + probeName, "rb") as fp:   #unPickling
        T1All.append(pickle.load(fp))
    print('finished reading data files...')



data=ProbeVars()
data.D = DAll.data
## correct the first D datapoint in time which is very high, dont know why!
data.D[:,0]=0

data.filterDA = filterDAAll.data
## correct the first D datapoint in time which is very high, dont know why!
data.filterDA[:,0]=0

## Temperature
data.T1= T1All.data
data.T1mean = np.mean(data.T1, axis=1)
data.T1prime = data.T1-data.T1mean[:,None]


## Velocity
data.U = UtnAll.data
data.u = data.U[:,:,0]
data.v = data.U[:,:,1]
data.w = data.U[:,:,2]

data.umean= np.mean(data.u, axis=1)
data.vmean= np.mean(data.v, axis=1)
data.wmean= np.mean(data.w, axis=1)

data.uprime = data.u-data.umean[:,None]
data.vprime = data.v-data.vmean[:,None]
data.wprime = data.w-data.wmean[:,None]

data.loctn = DAll.loctn
data.time= DAll.time

indexedProbes= ProbeIndex(DAll.loctn,[1e-3, 1e-9,0])
## choice of detector function


# data.detector = (np.square(ddt(np.multiply(data.uprime,data.vprime),data.time)))
# temp=highpass_filter(data.w,1/dt,1/(dt*50),2/(dt*50),1001)
# data.detector = np.square(ddt(temp,data.time))
# data.detector =(np.square(ddt(data.v,data.time)))
# data.detector =(np.square(ddt(data.w,data.time))+np.square(ddt(data.v,data.time)))
# data.detector =np.square( np.absolute(ddt(data.w,data.time)) + np.absolute(ddt(data.v,data.time)) )
# data.detector =(np.square(ddt(data.w,data.time))) # 
# data.detector = data.D/(data.umean[:,None]) 
# data.detector =( np.absolute(data.w) + np.absolute(data.v) )/(data.umean[:,None]) 
# data.detector =np.square( np.absolute(ddt(data.wprime,data.time)) + np.absolute(ddt(data.vprime,data.time)) + np.absolute(ddt(data.uprime,data.time)) )/(data.umean[:,None]) 

#data.detector = np.sqrt(data.filterDA)#/(data.umean[:,None]) #(np.sqrt(np.mean(np.square(data.wprime[:,None])))+np.sqrt(np.mean(np.square(data.vprime[:,None]))))
#data.detector = (np.square(ddt(np.multiply(data.uprime,data.vprime),data.time)))
data.detector=np.zeros(data.filterDA.shape)
## scaling detecors with freestream velocity
for j in np.arange(0,90,1):
    selectedProbesCurrentS = selectProbeIndexes(indexedProbes, [j, j], [0,200])
    selectedProbesCurrentSFreeStream = selectProbeIndexes(indexedProbes, [j, j], [20,20])
    
    
    data.detector[selectedProbesCurrentS,:]=np.sqrt(data.filterDA[selectedProbesCurrentS,:])/np.mean(data.umean[selectedProbesCurrentSFreeStream])
    #data.detector[selectedProbesCurrentS,:]= ((ddt(np.multiply(data.uprime[selectedProbesCurrentS,:],data.vprime[selectedProbesCurrentS,:]),data.time)))/np.mean(data.umean[selectedProbesCurrentSFreeStream])
    #data.detector[selectedProbesCurrentS,:] =np.absolute(ddt(data.wprime[selectedProbesCurrentS,:]+data.vprime[selectedProbesCurrentS,:]+data.uprime[selectedProbesCurrentS,:],data.time)) /np.mean(data.umean[selectedProbesCurrentSFreeStream])
# print(data.w.shape)
# print(data.detector.shape)
# print(data.D.shape)
# param1=np.absolute(ddt(data.w,data.time)) + np.absolute(ddt(data.v,data.time))/data.umean[:,None] 
# param2=np.square(param1)
# param3=np.square(param2)
# param4=np.square(param3)
# param52d=np.zeros(param2.shape)
# N=200
# for i in range(param52d.shape[0]):
#     #dFiltered2d[i,:]=np.convolve(data.detector[i,:], np.ones((N,))/N, mode='same')
#     param52d[i,:]=signal.convolve(param2[i,:], np.ones((N,))/N, mode='same',method='auto') # this is faster
# param5 = np.reshape(param52d,-1)
# pl.figure(1) 
# pl.hist(np.reshape(param1,-1)/np.max(np.reshape(param1,-1)),500,histtype='step',color='y',label="$param1$")
# pl.hist(np.reshape(param2,-1)/np.max(np.reshape(param2,-1)),500,histtype='step',color='r',label="$param2$")
# pl.hist(np.reshape(param3,-1)/np.max(np.reshape(param3,-1)),500,histtype='step',color='b',label="$param3$")
# pl.hist(np.reshape(param4,-1)/np.max(np.reshape(param4,-1)),500,histtype='step',color='k',label="$param4$")
# pl.hist(np.reshape(param5,-1)/np.max(np.reshape(param5,-1)),500,histtype='step',color='c',label="$param5$")
# pl.axis([0,1,0,1e5])
# legend = pl.legend(fontsize=14,numpoints=1,loc='upper right') 
# frame = legend.get_frame()
# frame.set_edgecolor('black')
# pl.draw()
## smoothing the signal
dx=0.01
dt=1e-5

N=200

C = 0.01

#################################### SELECTING PROBES

smin=0.25
smax=1
sIndexes=np.argwhere((data.loctn[:,0]>smin) & (data.loctn[:,0]<smax)  )


dFiltered2d=np.zeros((len(sIndexes),len(data.time)))
dFiltered2d2=np.zeros((len(sIndexes),len(data.time)))
dFiltered2d3=np.zeros((len(sIndexes),len(data.time)))
dFiltered2d4=np.zeros((len(sIndexes),len(data.time)))
for i in np.arange(sIndexes.shape[0]):
    iindex=sIndexes[i][0]
    #N=int(dx/data.umean[i]/dt)
    #print(N)
    #dFiltered2d[i,:]=np.convolve(data.detector[i,:], np.ones((N,))/N, mode='same')
    dFiltered2d[i,:]=signal.convolve(data.detector[iindex,:], np.ones((N,))/N, mode='same',method='auto') # this is faster
    dFiltered2d2[i,:]=signal.convolve(dFiltered2d[i,:], np.ones((N,))/N, mode='same',method='auto') # this is faster
    dFiltered2d3[i,:]=signal.convolve(dFiltered2d2[i,:], np.ones((N,))/N, mode='same',method='auto') # this is faster
    #dFiltered2d4[i,:]=generic_filter(data.detector[i,:], np.std, size=N)
    #dFiltered2d3[i,:]=signal.convolve(np.square(data.detector[i,:]), np.ones((N,))/N, mode='same',method='auto') # this is faster
    #dFiltered2d4[i,:]=np.sqrt(dFiltered2d3[i,:]-np.square(dFiltered2d[i,:]))
dFiltered = np.reshape(dFiltered2d,-1)

globalMax = np.amax(dFiltered)
ostuGlobalThresh = filters.threshold_otsu(dFiltered)
print(f'ostu thresh= {ostuGlobalThresh:5e}')
print('ostu thresh / max value = ', ostuGlobalThresh/globalMax )
print('global max=', globalMax)
print('C*global max=', C*globalMax)
# pl.figure(11)
# pl.plot(dFiltered)
# pl.show()

pl.figure(2) 
n1,x1,_=pl.hist(np.reshape(data.detector[data.detector<globalMax],-1),5000,histtype='step',density='true',color='y',label="$D$")#,cumulative=-1)
n2,x2,_=pl.hist(dFiltered,5000,histtype='step',color='k',density='true',label=r"$\displaystyle\widetilde{D}$")#,cumulative=-1)
print(n2)
print(np.diff(n1,2))
pl.figure(11)
pl.plot(np.diff(n1,2))
#pl.show()
print(x1[np.argmax(np.diff(n1,2))])
# pl.hist(data.detector[:,100],500,histtype='step',density='true',color='y',label="$D$")
# pl.hist(dFiltered2d[:,100],100,histtype='step',color='k',density='true',label=r"$\displaystyle\widetilde{D}$")
pl.gcf().subplots_adjust(bottom=0.15,left=0.15)
# pl.axis([0,1e9,0,5e-10])
legend = pl.legend(fontsize=14,numpoints=1,loc='upper right') 
frame = legend.get_frame()
frame.set_edgecolor('black')

# print(np.size(dFiltered[dFiltered>2e9]))
# print(np.size(dFiltered))
# density = stats.gaussian_kde(dFiltered[dFiltered>1e9])
# pl.plot(x,density(x))
pl.xlabel("S")
pl.ylabel("f(S)")
# pl.show()



selectedProbes = selectProbeIndexes(indexedProbes, [0, 200], [0,200 ])
selectedProbesFreeStream = selectProbeIndexes(indexedProbes, [0, 200], [20,20 ])
AvgDetectorFS=np.mean(data.detector[selectedProbesFreeStream,:])
selectedProbesTrail = selectProbeIndexes(indexedProbes, [89, 89], [0,200 ])
# sindex=36, s=0.299
ostuWThreshWallindex = np.zeros(21)
manualThreshWallindex = np.zeros(21)
for j in np.arange(0,21,1):
    print("wall index=",j)
    selectedProbesCurrentWallDistanceAll = selectProbeIndexes(indexedProbes, [0, 200], [j,j])
    selectedProbesCurrentWallDistanceLam = selectProbeIndexes(indexedProbes, [0, 36], [j,j])
    selectedProbesCurrentWallDistanceTurb = selectProbeIndexes(indexedProbes, [37, 200], [j,j])
    ostuWThreshWallindex[j] = filters.threshold_otsu(filterSignal(data.detector[selectedProbesCurrentWallDistanceAll,1:].flatten(),200,1))
    pl.subplot(211)
    maxAll=np.max(filterSignal(data.detector[selectedProbesCurrentWallDistanceAll,:].flatten(),200,1))
    maxLam=np.max(filterSignal(data.detector[selectedProbesCurrentWallDistanceLam,:].flatten(),200,1))
    maxTurb=np.max(filterSignal(data.detector[selectedProbesCurrentWallDistanceTurb,:].flatten(),200,1))
    maxBinvalue = np.min([maxAll,maxLam,maxTurb])
    countsAll,binsAll=np.histogram(filterSignal(data.detector[selectedProbesCurrentWallDistanceAll,:].flatten(),200,1),bins=1000,range=(0.0,maxBinvalue))
    countsLam,binsLam=np.histogram(filterSignal(data.detector[selectedProbesCurrentWallDistanceLam,:].flatten(),200,1),bins=1000,range=(0.0,maxBinvalue))
    countsTurb,binsTurb=np.histogram(filterSignal(data.detector[selectedProbesCurrentWallDistanceTurb,:].flatten(),200,1),bins=1000,range=(0.0,maxBinvalue))
    pl.hist(binsAll[:-1], binsAll, weights=countsAll/np.sum(countsAll),histtype='step',color='k',label=r"$\displaystyle\widetilde{D}$")
    pl.hist(binsLam[:-1], binsLam, weights=countsLam/np.sum(countsAll),histtype='step',color='b',label=r"$\displaystyle\widetilde{D}$")
    pl.hist(binsTurb[:-1], binsTurb, weights=countsTurb/np.sum(countsAll),histtype='step',color='r',label=r"$\displaystyle\widetilde{D}$")
    pl.axvline(ostuWThreshWallindex[j], color='k', linestyle='dashed', linewidth=1)
    pl.show()
    # # # pl.subplot(211)
    # # # pl.hist(filterSignal(data.detector[selectedProbesCurrentWallDistanceAll,:].flatten(),200,1),1000,histtype='step'
    # # # ,color='k',label=r"$\displaystyle\widetilde{D}$")#,cumulative=-1) #,density='true'
    # # # pl.hist(filterSignal(data.detector[selectedProbesCurrentWallDistanceLam,:].flatten(),200,1),1000,histtype='step'
    # # # ,color='b',label=r"$\displaystyle\widetilde{D}$")#,cumulative=-1)
    # # # pl.hist(filterSignal(data.detector[selectedProbesCurrentWallDistanceTurb,:].flatten(),200,1),1000,histtype='step'
    # # # ,color='r',label=r"$\displaystyle\widetilde{D}$")#,cumulative=-1)
    # # # pl.axvline(ostuWThreshWallindex[j], color='k', linestyle='dashed', linewidth=1)

    # # # pl.subplot(212)
    # # # pl.hist(filterSignal(data.detector[selectedProbesCurrentWallDistanceAll,:].flatten(),200,1),1000,histtype='step'
    # # # ,color='k',label=r"$\displaystyle\widetilde{D}$",cumulative=-1)
    # # # pl.hist(filterSignal(data.detector[selectedProbesCurrentWallDistanceLam,:].flatten(),200,1),1000,histtype='step'
    # # # ,color='b',label=r"$\displaystyle\widetilde{D}$",cumulative=-1)
    # # # pl.hist(filterSignal(data.detector[selectedProbesCurrentWallDistanceTurb,:].flatten(),200,1),1000,histtype='step'
    # # # ,color='r',label=r"$\displaystyle\widetilde{D}$",cumulative=-1)
    # # # #pl.axis([0,30000,0,1])
    # # #pl.show()
    manualThreshWallindex[j]=eval(input("Please enter the thresh value:"))
nprobes=len(selectedProbes)
indexedProbesSoretedByProbeNo=indexedProbes[np.argsort(indexedProbes[:,5])]
for i in range(nprobes):
    probe=selectedProbes[i]
    loctn=data.loctn[probe]
    ##print('probe location=',loctn)
    ## calculating detector at current freestream
    
    #indexprobeThisProbe=indexedProbes[indexedProbes[:,5]==probe].flatten()
    currentSIndex=indexedProbesSoretedByProbeNo[probe,0]#indexprobeThisProbe[0]
    currentWIndex=indexedProbesSoretedByProbeNo[probe,1]#indexprobeThisProbe[1]
    #selectedProbesCurrentFreeStream = selectProbeIndexes(indexedProbes, [currentSIndex, currentSIndex], [20,20 ])
    #selectedProbesCurrentWallDistance = selectProbeIndexes(indexedProbes, [0, 200], [currentWIndex,currentWIndex])
    #selectedProbesCurrentWS = selectProbeIndexes(indexedProbes, [currentSIndex, currentSIndex], [currentWIndex,currentWIndex])
    #RefDetectorFS=np.mean(data.detector[selectedProbesCurrentFreeStream,1:])
    #ostuWThresh = filters.threshold_otsu(filterSignal(data.detector[selectedProbesCurrentWallDistance,1:].flatten(),200,1))
    ##print('wall index=',currentWIndex)
    # # print('otsu thresh=',ostuWThreshWallindex[int(currentWIndex)])
    # # pl.subplot(211)
    # # pl.hist(filterSignal(data.detector[selectedProbesCurrentWS,:].flatten(),200,1),1000,histtype='step'
    # # ,color='k',density='true',label=r"$\displaystyle\widetilde{D}$")#,cumulative=-1)
    # # pl.axvline(ostuWThresh, color='k', linestyle='dashed', linewidth=1)
    # # pl.subplot(212)
    # # pl.hist(filterSignal(data.detector[selectedProbesCurrentWS,:].flatten(),200,1),1000,histtype='step'
    # # ,color='k',density='true',label=r"$\displaystyle\widetilde{D}$",cumulative=-1)
    # # #pl.axis([0,30000,0,1])
    # # pl.show()
    #print('Location=',loctn)
    time=data.time
    #N=int(dx/data.umean[probe]/dt)
    lamtu=lt(data.detector[probe,:],'thresh',manualThreshWallindex[int(currentWIndex)],refLam=data.detector[1,:]) #threshglob,ostuGlobalThresh,otsuloc
    gamma = np.append(gamma, [[loctn[0], loctn[1], loctn[2], np.mean(lamtu)]], axis = 0)

    # xlim=[np.min(time), np.max(time)]#[0.9, 1]#
    # ylim=[-2,2]
    # pl.figure
    
    # pl.title("s/L= " + str((loctn[0]/L).round(4)) + "d/L= " + str((loctn[1].round(5)/L).round(6)))
    # #pl.plot(time,Wprime,',',color='0.5')
    # pl.plot(time[np.where(lamtu==0)],data.u[probe,:][np.where(lamtu==0)],',',color='0.5')
    # pl.plot(time[np.where(lamtu==1)],data.u[probe,:][np.where(lamtu==1)],',r')
    # pl.plot(xlim,[np.mean(data.u[probe,:]),np.mean(data.u[probe,:])],'-.g',linewidth=2.0)
    # pl.plot(xlim,[np.mean(data.u[probe,:][np.where(lamtu==1)]),np.mean(data.u[probe,:][np.where(lamtu==1)])],'r',linewidth=2.0)
    # pl.plot(xlim,[np.mean(data.u[probe,:][np.where(lamtu==0)]),np.mean(data.u[probe,:][np.where(lamtu==0)])],color='0.5',linewidth=2.0)
    # # pl.plot(time[np.where(lamtu==0)],data.T1prime[probe,:][np.where(lamtu==0)],',',color='0.5')
    # # pl.plot(time[np.where(lamtu==1)],data.T1prime[probe,:][np.where(lamtu==1)],',r')
    # pl.xlabel('$t (s)$')
    # pl.ylabel('$u_t (m/s)$')
    # #pl.ylim(ylim)
    # pl.xlim(xlim)
    # pl.show()


    # xlim=[np.min(time), np.max(time)]#[0.9, 1]#
    # ylim=[-2,2]
    # pl.figure(figsize=(16,10))
    # pl.subplot(421)
    # pl.title("s/L= " + str((loctn[0]/L).round(4)) + "d/L= " + str((loctn[1].round(5)/L).round(6)))
    # #pl.plot(time,Wprime,',',color='0.5')
    # pl.plot(time[np.where(lamtu==0)],data.wprime[probe,:][np.where(lamtu==0)],',',color='0.5')
    # pl.plot(time[np.where(lamtu==1)],data.wprime[probe,:][np.where(lamtu==1)],',r')
    # # pl.plot(time[np.where(lamtu==0)],data.T1prime[probe,:][np.where(lamtu==0)],',',color='0.5')
    # # pl.plot(time[np.where(lamtu==1)],data.T1prime[probe,:][np.where(lamtu==1)],',r')

    # pl.ylim(ylim)
    # pl.xlim(xlim)

    # pl.subplot(423)
    # #pl.plot(time,Vprime,',',color='0.5')
    # pl.plot(time[np.where(lamtu==0)],data.vprime[probe,:][np.where(lamtu==0)],',',color='0.5')
    # pl.plot(time[np.where(lamtu==1)],data.vprime[probe,:][np.where(lamtu==1)],',r')
    # pl.ylim([-.4,0.4])
    # pl.xlim(xlim)

    # pl.subplot(425)
    # #pl.plot(time,Uprime,',',color='0.5')
    # pl.plot(time[np.where(lamtu==0)],data.uprime[probe,:][np.where(lamtu==0)],',',color='0.5')
    # pl.plot(time[np.where(lamtu==1)],data.uprime[probe,:][np.where(lamtu==1)],',r')
    # pl.ylim(ylim)
    # pl.xlim(xlim)
    # pl.subplot(427)
    # pl.plot(time[np.where(lamtu==0)],np.multiply(data.uprime,data.vprime)[probe,:][np.where(lamtu==0)],',',color='0.5')
    # pl.plot(time[np.where(lamtu==1)],np.multiply(data.uprime,data.vprime)[probe,:][np.where(lamtu==1)],',r')
    # pl.ylim([-.3,0.3])
    # pl.xlim(xlim)

    # pl.subplot(222)
    # pl.plot(data.uprime[probe,:][np.where(lamtu==0)],data.vprime[probe,:][np.where(lamtu==0)],',',color='0.8')
    # pl.plot(data.uprime[probe,:][np.where(lamtu==1)],data.vprime[probe,:][np.where(lamtu==1)],',r')
    
    # pl.subplot(224)
    # pl.plot(data.wprime[probe,:][np.where(lamtu==0)],data.vprime[probe,:][np.where(lamtu==0)],',',color='0.8')
    # pl.plot(data.wprime[probe,:][np.where(lamtu==1)],data.vprime[probe,:][np.where(lamtu==1)],',r')
    # pl.show()
GammaList, gammaarray=clusterProbes(gamma,[1e-3, 1e-9,0])


##print(gammaarray)
pl.figure(12)
pl.tricontour(gammaarray[:,0]/L, np.divide(gammaarray[:,1],gammaarray[:,3]), gammaarray[:,2], levels=np.arange(0.1,1,0.1), linewidths=1.0, colors='k' )
cntr = pl.tricontourf(gammaarray[:,0]/L, np.divide(gammaarray[:,1],gammaarray[:,3]), gammaarray[:,2], levels=50, linewidths=0.5, cmap="Reds")
pl.colorbar(cntr)
pl.xlabel("$s/L$")
pl.ylabel("$d/\delta_{99}$")

pl.figure(13)
pl.tricontour(gammaarray[:,0]/L, gammaarray[:,1], gammaarray[:,2], levels=np.arange(0.1,1,0.1), linewidths=1.0, colors='k')
cntr = pl.tricontourf(gammaarray[:,0]/L, gammaarray[:,1], gammaarray[:,2], levels=50, linewidths=0.5, cmap="Reds")
pl.colorbar(cntr)
pl.xlabel("$s/L$")
pl.ylabel("$d$")

pl.figure(14)
pl.tricontour(gammaarray[:,0]/L, gammaarray[:,1]/L, gammaarray[:,2], levels=np.arange(0.1,1,0.1), linewidths=1.0, colors='k')
cntr = pl.tricontourf(gammaarray[:,0]/L, gammaarray[:,1]/L, gammaarray[:,2], levels=50, linewidths=0.5, cmap="Reds")
pl.colorbar(cntr)
pl.xlabel("$s/L$")
pl.ylabel("$d/L$")
# pl.show()
pl.figure(7)
peakGamma=[]
for i in  range(0,len(GammaList),1): # [40,45,50,55,60,65,70,75,80,85,89]: #
    
    g = GammaList[i].gamma
    w = GammaList[i].wallDist/GammaList[i].delta99
    peakGamma.append([GammaList[i].s/L,np.max(g)])
    #print('i=',i)
    #print(str((GammaList[i].s/L).round(3)))  
    gsmooth = interp1d(w, g, kind='cubic',fill_value='extrapolate')
    wmid=w[:-1] + np.diff(w) / 2
    wnew = np.sort(np.concatenate((w,wmid)))
    
    gnew = gsmooth(wnew)
    pl.plot(w,g,label="$s/L=$"+str((GammaList[i].s/L).round(2)))
    legend = pl.legend(fontsize=14,numpoints=1,loc='upper right') 
    frame = legend.get_frame()
    frame.set_edgecolor('black')
    pl.xlabel("$d/\delta_{99}$")
    pl.ylabel("$\gamma$")
    pl.xlim([0,1.5])
    pl.ylim([0, 1])

pl.figure(8)
peakGamma=np.asarray(peakGamma)
##print(peakGamma)
pl.plot(peakGamma[:,0],peakGamma[:,1],'k',linewidth=3.0,label='current study')
popt, pcov=curve_fit(general_gamma, peakGamma[:,0], peakGamma[:,1],p0=[0.64,1.5])
##print(popt)
pl.plot(np.arange(0.3,2,0.001),general_gamma(np.arange(0.3,2,0.001),*popt),'--r',linewidth=2,label='$1-exp(-5\eta^3)$')
legend = pl.legend(fontsize=14,numpoints=1,loc='upper left') 
frame = legend.get_frame()
frame.set_edgecolor('black')
pl.xlabel("$s/L$")
pl.ylabel("$\gamma_{max}$")


pl.figure(9)
##plot selected lines
for i in  [40,45,50,55,60,65,70,75,80,85,89]: # range(0,len(GammaList),1): # 
    
    g = GammaList[i].gamma
    w = np.asarray(GammaList[i].wallDist)/GammaList[i].delta99
    
    #print('i=',i)
    #print(str((GammaList[i].s/L).round(3)))  
    gsmooth = interp1d(w, g, kind='cubic',fill_value='extrapolate')
    wmid=w[:-1] + np.diff(w) / 2
    wnew = np.sort(np.concatenate((w,wmid)))
    
    gnew = gsmooth(wnew)
    pl.plot(w,g,label="$s/L=$"+str((GammaList[i].s/L).round(2)))
    legend = pl.legend(fontsize=14,numpoints=1,loc='upper right') 
    frame = legend.get_frame()
    frame.set_edgecolor('black')
    pl.xlabel("$d/\delta_{99}$")
    pl.ylabel("$\gamma$")
    pl.xlim([0,1.5])
    pl.ylim([0, 1])

pl.figure(20)
##plot selected lines
for i in  [40,45,50,55,60,65,70,75,80,85,89]: # range(0,len(GammaList),1): # 
    
    g = GammaList[i].gamma
    w = np.asarray(GammaList[i].wallDist)/L
    
    #print('i=',i)
    #print(str((GammaList[i].s/L).round(3)))  
    gsmooth = interp1d(w, g, kind='cubic',fill_value='extrapolate')
    wmid=w[:-1] + np.diff(w) / 2
    wnew = np.sort(np.concatenate((w,wmid)))
    
    gnew = gsmooth(wnew)
    pl.plot(w,g,label="$s/L=$"+str((GammaList[i].s/L).round(2)))
    legend = pl.legend(fontsize=14,numpoints=1,loc='upper right') 
    frame = legend.get_frame()
    frame.set_edgecolor('black')
    pl.xlabel("$d/L$")
    pl.ylabel("$\gamma$")
    pl.xlim([0,0.03])
    pl.ylim([0, 1])


pl.show()

if yes_or_no('Would you like to save Gamma data?'):
    filename="GammaList_"+date_time.strftime("%b-%d-%Y-%H_%M_%S")
    print(filename)
    with open(filename, "wb") as fp:   #Pickling
        pickle.dump(GammaList, fp)




#selectedProbes = selectProbes(DAll.loctn, [0.34, 0.7], [0, 1e-4])

#print(test.U.shape)

##### MEthods:
##### 1-Otsu
##### 2-Adpative
##### 3-Aribitrary fraction
method=2
for i in range(nprobes):
    probe=selectedProbes[i]
    print("location= ", DAll.getProbe(probe).loctn.round(8))
    loctn=DAll.getProbe(probe).loctn
    time=DAll.getProbe(probe).time
    D=DAll.getProbe(probe).data #
    U=UtnAll.getProbe(probe).data[:,0]
    V=UtnAll.getProbe(probe).data[:,1]
    W=UtnAll.getProbe(probe).data[:,2]
    Umean=np.mean(U)
    Vmean=np.mean(V)
    Wmean=np.mean(W)
    Uprime = U - Umean
    Vprime = V - Vmean
    Wprime = W - Wmean
    #D=(np.square(ddt(W,time))) # 
    #D= (np.square(ddt(np.multiply(U,V),time)))
    #D=Allprobes.getProbe(probes).filterDA
    ##Dtemp=highpass_filter(W,1/dt,1/(dt*50),2/(dt*50),1001)
    ##D=np.square(ddt(Dtemp,time))
    detector=D
    ##smoothing the signal
    
    gamma = np.append(gamma, [[loctn[0], loctn[1], loctn[2], np.mean(lt)]], axis = 0)
    xlim=[np.min(time), np.max(time)]#[0.9, 1]#
    ylim=[-2,2]

    # pl.figure(figsize=(8,10))
    # pl.subplot(411)
    # pl.title("s/L= " + str((loctn[0]/L).round(4)) + "d/L= " + str((loctn[1].round(5)/L).round(6)))
    # #pl.plot(time,Wprime,',',color='0.5')
    # pl.plot(time[np.where(lt==0)],Wprime[np.where(lt==0)],',',color='0.5')
    # pl.plot(time[np.where(lt==1)],Wprime[np.where(lt==1)],',r')
    # pl.ylim(ylim)
    # pl.xlim(xlim)

    # pl.subplot(412)
    # #pl.plot(time,Vprime,',',color='0.5')
    # pl.plot(time[np.where(lt==0)],Vprime[np.where(lt==0)],',',color='0.5')
    # pl.plot(time[np.where(lt==1)],Vprime[np.where(lt==1)],',r')
    # pl.ylim([-.4,0.4])
    # pl.xlim(xlim)

    # pl.subplot(413)
    # #pl.plot(time,Uprime,',',color='0.5')
    # pl.plot(time[np.where(lt==0)],Uprime[np.where(lt==0)],',',color='0.5')
    # pl.plot(time[np.where(lt==1)],Uprime[np.where(lt==1)],',r')
    # pl.ylim(ylim)
    # pl.xlim(xlim)
    # pl.subplot(414)
    # pl.plot(Uprime[np.where(lt==0)],Vprime[np.where(lt==0)],',',color='0.8')
    # pl.plot(Uprime[np.where(lt==1)],Vprime[np.where(lt==1)],',r')
    
    # pl.show()
    
#print('##############################MAX')
#print(np.max(maxglobal))
GammaList=clusterProbes(gamma,[1e-3, 1e-9,0])
#with open("GammaList.txt", "wb") as fp:   #Pickling
#    pickle.dump(GammaList, fp)


pl.figure(7)
for i in range(0,len(GammaList),6):
    
    g = GammaList[i].gamma
    w = GammaList[i].wallDist
       
    gsmooth = interp1d(w, g, kind='cubic',fill_value='extrapolate')
    wmid=w[:-1] + np.diff(w) / 2
    wnew = np.sort(np.concatenate((w,wmid)))
    
    gnew = gsmooth(wnew)
    pl.plot(w,g,label="$s/L=$"+str((GammaList[i].s/L).round(2)))
    legend = pl.legend(fontsize=14,numpoints=1,loc='upper right') 
    frame = legend.get_frame()
    frame.set_edgecolor('black')
    pl.xlabel("Wall Distance")
    pl.ylabel("$\gamma$")
    pl.xlim([0,0.015])
    pl.ylim([0, 1])
pl.show()
alaki
    
pl.figure(1)
pl.plot(time,U-Umean,'-',label="$u^\prime_t$",color = '0.75')
pl.plot(time,V-Vmean,'--b',label="$u^\prime_n$")
pl.plot(time,W-Wmean,'-.g',label="$w^\prime$")
pl.xlabel("time, t")
pl.ylabel("velocity fluctuations")
pl.gcf().subplots_adjust(bottom=0.2,left=0.2)
legend = pl.legend(fontsize=14,numpoints=1,loc='upper right') 
frame = legend.get_frame()
frame.set_edgecolor('black')
pl.axis([0,np.amax(time)*1.1,-7,7])
###plotting detector
pl.figure(2)
pl.plot(time,D,'-y',label="$D$")
pl.plot(time[int(N/2):int(-N/2+1)],dFiltered,'k',label=r"$\displaystyle\widetilde{D}$")
pl.xlabel("time, t")
pl.ylabel("criterion function")
pl.gcf().subplots_adjust(bottom=0.2,left=0.2)
legend = pl.legend(fontsize=14,numpoints=1,loc='upper right') 
frame = legend.get_frame()
frame.set_edgecolor('black')

pl.figure(3) 
pl.hist(D,1000,histtype='step',color='y',label="$D$")
pl.hist(dFiltered,1000,histtype='step',color='k',label=r"$\displaystyle\widetilde{D}$")
pl.xlabel("S")
pl.ylabel("f(S)")
pl.gcf().subplots_adjust(bottom=0.15,left=0.15)
pl.axis([0,1e9,0,100])
legend = pl.legend(fontsize=14,numpoints=1,loc='upper right') 
frame = legend.get_frame()
frame.set_edgecolor('black')

####### plotting discrimination over velocity signal



fig4=pl.figure(figsize=(10,4))
pl.plot(time,U-Umean,'-',label="$u^\prime_t$",color = '0.75')
pl.plot(time,V-Vmean,'--b',label="$u^\prime_n$")
pl.plot(time,W-Wmean,'-.g',label="$w^\prime$")
pl.xlabel("time, t")
pl.ylabel("velocity fluctuations")
pl.gcf().subplots_adjust(bottom=0.2,left=0.2)
legend = pl.legend(fontsize=14,numpoints=1,loc='upper right') 
frame = legend.get_frame()
frame.set_edgecolor('black')
pl.axis([0,np.amax(time),-7,7])

pl.plot(time[int(N/2):int(-N/2+1)],lt,'-r',label="LT")
fig4.tight_layout()
####################################
# ###plotting detector
####################################
fig5=pl.figure(figsize=(10,4))
#pl.figure(5)
pl.plot(time,D,'-y',label="$D$")
pl.plot(time[int(N/2):int(-N/2+1)],dFiltered,'k',label=r"$\displaystyle\widetilde{D}$")
pl.xlabel("time, t")
pl.ylabel("criterion function")
pl.gcf().subplots_adjust(bottom=0.2,left=0.2)
legend = pl.legend(fontsize=14,numpoints=1,loc='upper right') 
frame = legend.get_frame()
frame.set_edgecolor('black')
pl.plot(time[int(N/2):int(-N/2+1)],(lt/5)*np.amax(D)/2,'-r',label="LT")
pl.axis([0,np.amax(time),0,np.amax(D)])
#pl.plot(time[int(N/2):int(-N/2+1)],lt_adaptive[0,:]*np.amax(dFiltered)*2/11,'-k',label="LT") 
fig5.tight_layout()

############################
##filtered signal
############################
pl.figure(6)
pl.plot(time,W-Wmean,'-y',label="$W-Wmean$")
pl.plot(time,highpass_filter(W-Wmean,1/dt,1/(dt*50),2/(dt*50),1001),'k',label=r"$\displaystyle\widetilde{D}$")

#############################
### Energy Spectra
#############################
Fs=1/dt
Nw=8192
f , pAvg = signal.welch(W, Fs,scaling='density',noverlap=0,nperseg=Nw)

pl.figure(7)
pl.gcf().set_size_inches(7.2, 4.4, forward=True)
#pl.plot(kEta, p1Avg, 'k')
pl.plot(f, pAvg, 'k')
#pl.axis([0.1,50,0.01,1e6])
#pl.plot([2 , 20],[130000, 10**(np.log10(130000)-5/3)],'--k')


pl.show()
# print()
alaki
cc=1
while cc:
    cc=eval(input("do you want to continue:"))
    
    if cc>0 :
        probes=eval(input("Please enter the probe number:"))
        #print("location= ", Allprobes.getProbe(probes).loctn.round(10))
        pl.figure(1)
        pl.plot(Allprobes.getProbe(probes).time,Allprobes.getProbe(probes).Utn[:,1],'-k',label="Ut_y")
        pl.plot(Allprobes.getProbe(probes).time,Allprobes.getProbe(probes).Utn[:,2],'-.k',label="Ut_w")
        pl.plot(Allprobes.getProbe(probes).time,Allprobes.getProbe(probes).Utn[:,0],'-b',label="Ut_x")
        legend = pl.legend(fontsize=20,numpoints=1) 
        frame = legend.get_frame()
        frame.set_edgecolor('black')
        pl.xlabel("time, t")
        pl.ylabel("Utn")
        pl.gcf().subplots_adjust(bottom=0.15,left=0.15)
        pl.figure(2)
        pl.plot(Allprobes.getProbe(probes).time,Allprobes.getProbe(probes).D,'-k',label="D")
        legend = pl.legend(fontsize=20,numpoints=1) 
        frame = legend.get_frame()
        frame.set_edgecolor('black')
        pl.xlabel("time, t")
        pl.ylabel("D")
        pl.gcf().subplots_adjust(bottom=0.15,left=0.15)
        pl.show()

ifsave=eval(input("do you want to save this probe:"))
if ifsave:
    fout = open("F:/PostProcess/ccc4/postProcessing/test/0/Utn",'a')        
    data=np.c_[Allprobes.getProbe(probes).time,Allprobes.getProbe(probes).Utn[:,0],Allprobes.getProbe(probes).Utn[:,1],Allprobes.getProbe(probes).Utn[:,2]]

    #np.savetxt(fout,St,fmt="%10.5f ")
    np.savetxt(fout,data,header=" Probe "+str(probes) + " " + str(Allprobes.getProbe(probes).loc),fmt="%15.12f %15.12f %15.12f %15.12f")
    fout.close()
    fout = open("F:/PostProcess/ccc4/postProcessing/test/1e-05/D",'a')        
    data=np.c_[Allprobes.getProbe(probes).time,Allprobes.getProbe(probes).D]

    #np.savetxt(fout,St,fmt="%10.5f ")
    np.savetxt(fout,data,header=" Probe "+str(probes) + " " + str(Allprobes.getProbe(probes).loc),fmt="%15.12f %15.12f")
    fout.close()
    
    fout = open("F:/PostProcess/ccc4/postProcessing/test/1e-05/filterDA",'a')        
    data=np.c_[Allprobes.getProbe(probes).time,Allprobes.getProbe(probes).filterDA]

    #np.savetxt(fout,St,fmt="%10.5f ")
    np.savetxt(fout,data,header=" Probe "+str(probes) + " " + str(Allprobes.getProbe(probes).loc),fmt="%15.12f %15.12f")
    fout.close()
    
    fout = open("F:/PostProcess/ccc4/postProcessing/test/1e-05/filterDG",'a')        
    data=np.c_[Allprobes.getProbe(probes).time,Allprobes.getProbe(probes).filterDG]

    #np.savetxt(fout,St,fmt="%10.5f ")
    np.savetxt(fout,data,header=" Probe "+str(probes) + " " + str(Allprobes.getProbe(probes).loc),fmt="%15.12f %15.12f")
    fout.close()
alaki


#