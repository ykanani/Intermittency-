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
from scipy.interpolate import spline
#os.chdir("D:\Dropbox\PhD\Python")
import easygui
from matplotlib import rc
from probeAnalysis import Probes
from probeAnalysis import Gamma
from probeAnalysis import selectProbes
from probeAnalysis import processFile
from probeAnalysis import clusterProbes
from probeAnalysis import highpass_filter

#####GOOD FONTS!!!! remove comment
rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':20})
rc('text', usetex=True)
os.chdir(os.path.dirname(os.path.realpath(__file__)))
print(os.getcwd())

os.chdir('F:/PostProcess/ccc2')
os.listdir()
print(os.getcwd())

##input parameters
nu = 1.55e-5
#nu =7.7873794212218649517684887459807e-6 #for scaled length and velocity
C=0.4976 #chord
Uin=4
T=C/Uin
dt=1e-5

#####################
##################### Main
#####################


ndelta=21  # number of probes at each station (wall normal)
ns=90      # number of stations
probeInfo=[0,ndelta*ns-1,1] # [first,last,interval]
zlocs=['-0.05', '-0.04']#, '-0.03', '-0.02', '-0.01', '0', '0.01', '0.02', '0.03', '0.04', '0.05']
gamma = np.empty([0,4]) # culumns: s,w,z,gamma
for k in range(len(zlocs)):
    probeName=  "p_ss_BL_z" + zlocs[k]# "test"  # 
    print('reading data file',"probes-" + probeName)
    with open("probes-" + probeName, "rb") as fp:   #Pickling
        Allprobes=pickle.load(fp)
    
    print(Allprobes.loctn[:,0].round(4))
    selectedProbes = selectProbes(Allprobes.loctn, [0, 0.7], [0, 100])
    nprobes=len(selectedProbes)
    #print(test.U.shape)
    
    ##### MEthods:
    ##### 1-Otsu
    ##### 2-Adpative
    ##### 3-Aribitrary fraction
    method=3
    for i in range(nprobes):
        probe=selectedProbes[i]
        print("location= ", Allprobes.getProbe(probe).loctn.round(10))
        loctn=Allprobes.getProbe(probe).loctn
        time=Allprobes.getProbe(probe).time
        U=Allprobes.getProbe(probe).Utn[:,0]
        V=Allprobes.getProbe(probe).Utn[:,1]
        W=Allprobes.getProbe(probe).Utn[:,2]
        Umean=np.mean(Allprobes.getProbe(probe).Utn[:,0])
        Vmean=np.mean(Allprobes.getProbe(probe).Utn[:,1])
        Wmean=np.mean(Allprobes.getProbe(probe).Utn[:,2])

        D=Allprobes.getProbe(probe).D #
        #D=(np.square(ddt(W,time))) # 
        #D= (np.square(ddt(np.multiply(U,V),time)))
        #D=Allprobes.getProbe(probes).filterDA
        ##Dtemp=highpass_filter(W,1/dt,1/(dt*50),2/(dt*50),1001)
        ##D=np.square(ddt(Dtemp,time))

        ##smoothing the signal
        N=100
        DSmooth=np.convolve(D, np.ones((N,))/N, mode='valid')
        ## Laminar/Turbulent Discrimination
        if method==1:
            ### Otsu
            thresh_otsu= filters.threshold_otsu(DSmooth)
            lt_otsu=np.where(DSmooth>thresh_otsu, 1, 0) # Laminar_Turbulent, Ostu
            lt = lt_otsu
            print("Otsu Threshhold={:.2e}".format(thresh_otsu))
            #pl.plot(time[int(N/2):int(-N/2+1)],lt_otsu,'-r',label="LT")
        elif method==2:
            # ### Adaptive
            block_size = 50001
            print(np.tile(DSmooth,(2,1)).shape)
            lt_adaptive = filters.threshold_adaptive(np.tile(DSmooth,(2,1)), block_size, offset=100)*1
            lt = lt_adaptive[0,:]
            print(lt_adaptive[0,:].shape)
            #pl.plot(time[int(N/2):int(-N/2+1)],lt,'-k',label="LT") 
        elif method==3:
            ### Arbitrary fraction of Max D
            C=0.01
            lt_cmax=np.where(DSmooth>C*np.max(DSmooth), 1, 0) # Laminar_Turbulent, C*max(D)
            lt = lt_cmax
            #pl.plot(time[int(N/2):int(-N/2+1)],lt_cmax,'-b',label="LT")
        gamma = np.append(gamma, [[loctn[0], loctn[1], loctn[2], np.mean(lt)]], axis = 0)
        #gamma[i,0]=loctn[0]
        #gamma[i,1]=loctn[1]
        #gamma[i,2]=loctn[2]
        #gamma[i,3]=np.mean(lt)
    
    
print("gamma")
print(gamma)  
print("arg sort s")

GammaList=clusterProbes(gamma,[1e-3, 1e-9,0])
with open("GammaList.txt", "wb") as fp:   #Pickling
    pickle.dump(GammaList, fp)


pl.figure(7)
for i in range(len(GammaList)):
    print(GammaList[i].wallDist)
    pl.plot(GammaList[i].wallDist,GammaList[i].gamma,label="$\gamma$")
legend = pl.legend(fontsize=14,numpoints=1,loc='upper right') 
frame = legend.get_frame()
frame.set_edgecolor('black')
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
pl.plot(time[int(N/2):int(-N/2+1)],DSmooth,'k',label=r"$\displaystyle\widetilde{D}$")
pl.xlabel("time, t")
pl.ylabel("criterion function")
pl.gcf().subplots_adjust(bottom=0.2,left=0.2)
legend = pl.legend(fontsize=14,numpoints=1,loc='upper right') 
frame = legend.get_frame()
frame.set_edgecolor('black')

pl.figure(3) 
pl.hist(D,1000,histtype='step',color='y',label="$D$")
pl.hist(DSmooth,1000,histtype='step',color='k',label=r"$\displaystyle\widetilde{D}$")
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
pl.plot(time[int(N/2):int(-N/2+1)],DSmooth,'k',label=r"$\displaystyle\widetilde{D}$")
pl.xlabel("time, t")
pl.ylabel("criterion function")
pl.gcf().subplots_adjust(bottom=0.2,left=0.2)
legend = pl.legend(fontsize=14,numpoints=1,loc='upper right') 
frame = legend.get_frame()
frame.set_edgecolor('black')
pl.plot(time[int(N/2):int(-N/2+1)],(lt/5)*np.amax(D)/2,'-r',label="LT")
pl.axis([0,np.amax(time),0,np.amax(D)])
#pl.plot(time[int(N/2):int(-N/2+1)],lt_adaptive[0,:]*np.amax(DSmooth)*2/11,'-k',label="LT") 
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
        print("location= ", Allprobes.getProbe(probes).loctn.round(10))
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