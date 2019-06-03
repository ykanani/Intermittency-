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




ndelta=21  # number of probes at each station (wall normal)
ns=90      # number of stations
probeInfo=[0,ndelta*ns-1,1] # [first,last,interval]
zlocs=['-0.05', '-0.04']#, '-0.03', '-0.02', '-0.01', '0', '0.01', '0.02', '0.03', '0.04', '0.05']
gamma = np.empty([0,4]) # culumns: s,w,z,gamma
for k in range(len(zlocs)):
    probeName=  "p_ss_BL_z" + zlocs[k]# "test"  # 
    
    ## load from file
    with open("probes-" + probeName, "rb") as fp:   #Pickling
        Allprobes = pickle.load(fp)
    #print(Allprobes.loc)
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




### read processed gamma        
with open("GammaList.txt", "rb") as fp:   # Unpickling
    GammaList = pickle.load(fp)

pl.figure(7)
for i in range(40,len(GammaList),5):
    print(GammaList[i].wallDist)
    pl.plot(GammaList[i].wallDist,GammaList[i].gamma,label="$s/L=$"+str((GammaList[i].s/C).round(2)))
legend = pl.legend(fontsize=14,numpoints=1,loc='upper right') 
frame = legend.get_frame()
frame.set_edgecolor('black')
pl.xlabel("Wall Distance")
pl.ylabel("$\gamma$")
pl.show()