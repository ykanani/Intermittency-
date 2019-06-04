
import csv
# reads the data from multiple folders
#import statistics
import math
import gc
import os
import pickle
import re
import sys
import easygui
import matplotlib.pyplot as pl
import numpy as np
from matplotlib import rc
from scipy import signal
from scipy.interpolate import spline
from skimage import filters

from probeAnalysis import (Gamma, Probes, clusterProbes, processFile,
                           selectProbes)

os.chdir('F:/PostProcess/ccc2')
os.listdir()

print(os.getcwd())
##input parameters
nu = 1.55e-5
test= 1
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
zlocs=['-0.05', '-0.04', '-0.03', '-0.02', '-0.01', '0', '0.01', '0.02', '0.03', '0.04', '0.05']
gamma = np.empty([0,4]) # culumns: s,w,z,gamma
for k in range(len(zlocs)):
    probeName=  "p_ss_BL_z" + zlocs[k]# "test"  #
    
    # D = processFile(probeInfo,probeName,initTime= "0.43001",varName="D",nComp=1)
    # with open("postProcessing/" + D.varName + "_" + probeName, "wb") as fp:   #Pickling
    #     pickle.dump(D, fp)
    # del D
    # gc.collect()

    # p = processFile(probeInfo,probeName,initTime= "0.43",varName="p",nComp=1)
    # with open("postProcessing/" + p.varName + "_" + probeName, "wb") as fp:   #Pickling
    #     pickle.dump(p, fp)
    # del p
    # gc.collect()

    # T1 = processFile(probeInfo, probeName,initTime= "0.43",varName="T1",nComp=1)
    # with open("postProcessing/" + T1.varName + "_" + probeName, "wb") as fp:   #Pickling
    #     pickle.dump(T1, fp)
    # del T1
    

    # T2 = processFile(probeInfo,probeName,initTime= "0.43",varName="T2",nComp=1)
    # with open("postProcessing/" + T2.varName + "_" + probeName, "wb") as fp:   #Pickling
    #     pickle.dump(T2, fp)
    # del T2
    # gc.collect()

    # U = processFile(probeInfo,probeName,initTime= "0.43",varName="U",nComp=3)
    # with open("postProcessing/" + U.varName + "_" + probeName, "wb") as fp:   #Pickling
    #     pickle.dump(U, fp)
    # del U
    # gc.collect()

    Utn = processFile(probeInfo,probeName,initTime= "0.43",varName="Utn",nComp=3)
    with open("postProcessing/" + Utn.varName + "_" + probeName, "wb") as fp:   #Pickling
        pickle.dump(Utn, fp)
    del Utn
    gc.collect()

    # filterDA = processFile(probeInfo,probeName,initTime= "0.43001",varName="filterDA",nComp=1)
    # with open("postProcessing/" + filterDA.varName + "_" + probeName, "wb") as fp:   #Pickling
    #     pickle.dump(filterDA, fp)
    # del filterDA
    # gc.collect()

    # filterDG = processFile(probeInfo,probeName,initTime= "0.43001",varName="filterDG",nComp=1)
    # with open("postProcessing/" + filterDG.varName + "_"  + probeName, "wb") as fp:   #Pickling
    #     pickle.dump(filterDG, fp)
    # del filterDG 
    # gc.collect()

    # filterQA = processFile(probeInfo,probeName,initTime= "0.43001",varName="filterQA",nComp=1)
    # with open("postProcessing/" + filterQA.varName + "_" + probeName, "wb") as fp:   #Pickling
    #     pickle.dump(filterQA, fp)
    # del filterQA
    # gc.collect()

    # gradU = processFile(probeInfo,probeName,initTime= "0.43001",varName="gradU",nComp=9)
    # with open("postProcessing/" + gradU.varName + "_" + probeName, "wb") as fp:   #Pickling
    #     pickle.dump(gradU, fp)
    # del gradU
    # gc.collect()

    # Q = processFile(probeInfo,probeName,initTime= "0.43001",varName="Q",nComp=1)
    # with open("postProcessing/" + Q.varName + "_" + probeName, "wb") as fp:   #Pickling
    #     pickle.dump(Q, fp)
    # del Q
    # gc.collect()

    # vorticity = processFile(probeInfo,probeName,initTime= "0.43001",varName="vorticity",nComp=1)
    # with open("postProcessing/" + vorticity.varName + "_" + probeName, "wb") as fp:   #Pickling
    #     pickle.dump(vorticity, fp)
    # gc.collect()

    

    
#
