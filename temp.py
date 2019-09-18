for i in range(nprobes):
    probe=selectedProbes[i]
    loctn=data.loctn[probe]
    print('probe location=',loctn)
    ## calculating detector at current freestream
    
    #indexprobeThisProbe=indexedProbes[indexedProbes[:,5]==probe].flatten()
    currentSIndex=indexedProbesSoretedByProbeNo[probe,0]#indexprobeThisProbe[0]
    currentWIndex=indexedProbesSoretedByProbeNo[probe,1]#indexprobeThisProbe[1]
    selectedProbesCurrentFreeStream = selectProbeIndexes(indexedProbes, [currentSIndex, currentSIndex], [20,20 ])
    selectedProbesCurrentWallDistance = selectProbeIndexes(indexedProbes, [0, 200], [currentWIndex,currentWIndex])
    selectedProbesCurrentWS = selectProbeIndexes(indexedProbes, [currentSIndex, currentSIndex], [currentWIndex,currentWIndex])
    RefDetectorFS=np.mean(data.detector[selectedProbesCurrentFreeStream,1:])
    #ostuWThresh = filters.threshold_otsu(filterSignal(data.detector[selectedProbesCurrentWallDistance,1:].flatten(),200,1))
    print('wall index=',currentWIndex)
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
    lamtu=lt(data.detector[probe,:],'thresh',manualThreshWallindex[int(currentWIndex)]) #threshglob,ostuGlobalThresh,otsuloc ,refLam=data.detector[1,:]
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


print(gammaarray)
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
print(peakGamma)
pl.plot(peakGamma[:,0],peakGamma[:,1],'k',linewidth=3.0,label='current study')
popt, pcov=curve_fit(general_gamma, peakGamma[:,0], peakGamma[:,1],p0=[0.64,1.5])
print(popt)
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

filename="GammaList_"+date_time.strftime("%b-%d-%Y-%H_%M_%S")
print(filename)
with open(filename, "wb") as fp:   #Pickling
    pickle.dump(GammaList, fp)