import math
import gc
import os
import pickle
import re
import sys
#import easygui
import matplotlib.pyplot as pl
import numpy as np

from probeAnalysis import (Gamma, ProbeVars, clusterProbes, processFile,processWallProbe,selectProbes)

os.chdir('I:/PostProcess/ccc9')
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

probeName = "p_ss_spanwise_d0"
processWallProbe(probeName,initTime= "0.35",varName="line1_U_Utn_vorticity.xy",nComp=1)



