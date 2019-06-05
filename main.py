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
#import easygui
from matplotlib import rc
from probeAnalysis import Probes
from probeAnalysis import Gamma
from probeAnalysis import selectProbes
from probeAnalysis import selectProbeIndexes
from probeAnalysis import processFile
from probeAnalysis import clusterProbes
from probeAnalysis import highpass_filter
from probeAnalysis import ProbeIndex
from scipy.interpolate import interp1d
import matplotlib as mpl; print(mpl.font_manager.get_cachedir())
#####GOOD FONTS!!!! remove comment
# rc('font', **{'family': 'serif', 'serif': ['Computer Modern'],'size':20})
# rc('text', usetex=True)
os.chdir(os.path.dirname(os.path.realpath(__file__)))
print(os.getcwd())

os.chdir('I:/PostProcess/ccc2')
os.listdir()
print(os.getcwd())

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
zlocs=['-0.05', '-0.04']#, '-0.03', '-0.02', '-0.01', '0', '0.01', '0.02', '0.03', '0.04', '0.05']
gamma = np.empty([0,4]) # culumns: s,w,z,gamma
for k in range(len(zlocs)):
    probeName=  "p_ss_BL_z" + zlocs[k]# "test"  # 
    print('reading data file',"D" + "_" + probeName)
    with open("postProcessing/" + "D" + "_" + probeName, "rb") as fp:   #unPickling
        DAll=pickle.load(fp)
    print('reading data file',"Utn" + "_" + probeName)
    with open("postProcessing/" + "Utn" + "_" + probeName, "rb") as fp:   #unPickling
        UtnAll=pickle.load(fp)
    print('finished reading data file',"probes-" + probeName)
    print(DAll.loctn[:,0].round(4))
    problocs=np.insert(DAll.loctn,[2],np.ones([len(DAll.loctn[:,0]),1]),axis=1)
    indexedProbes= ProbeIndex(DAll.loctn,[1e-3, 1e-9,0])
    
    selectedProbes = selectProbeIndexes(indexedProbes, [80, 80], [0,25 ])
    
    #selectedProbes = selectProbes(DAll.loctn, [0.34, 0.7], [0, 1e-4])
    nprobes=len(selectedProbes)
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
        N=100
        maxglobal=[]
        DSmooth=np.convolve(detector, np.ones((N,))/N, mode='valid')
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
            
            maxglobal.append(np.max(DSmooth)/np.square(Umean))            #lt_cmax=np.where(DSmooth>C*np.max(DSmooth), 1, 0) # Laminar_Turbulent, C*max(D)
            lt_cmax=np.where(DSmooth>C*np.max(DSmooth),1,0)#0.01*24473540.92065062,1,0)#>0.01*269112747.3345178,1,0) ##6419964748.39, 1, 0) # Laminar_Turbulent, C*max(D)
            lt = lt_cmax
            #pl.plot(time[int(N/2):int(-N/2+1)],lt_cmax,'-b',label="LT")
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
    pl.plot(wnew,gnew,label="$s/L=$"+str((GammaList[i].s/L).round(2)))
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