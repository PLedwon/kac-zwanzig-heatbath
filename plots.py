#!/usr/bin/python3
import numpy as np
import glob
import matplotlib.pyplot as plt
import math
import scipy
from scipy.optimize import curve_fit

#def moving_average(a, n=3) :
#    ret = np.cumsum(a, dtype=float)
#    ret[n:] = ret[n:] - ret[:-n]
#    return ret[n - 1:] / n


if not glob.glob('../data/*.npz'):

    resultList = glob.glob('/users/stud/ledwon/Documents/npzFiles/*.np[yz]')
    
    data=np.load(resultList[0])
    varQ=np.zeros(np.size(data['squaredQ']))
    varP=np.zeros(np.size(data['squaredP']))
    aveQ=np.zeros(np.size(data['aveQ']))
    aveP=np.zeros(np.size(data['aveP']))
    K=np.zeros(np.size(data['K']))
    timesteps=data['timesteps']
    t1=data['t1']
    dt=data['dt']
    gamma=data['gamma']
    maxEError=0
    errorFileCount=0
    
    
    errorbarCount = 100
    indexSkipValue = int(1.0/dt * t1/float(errorbarCount))
    timestepsErr=timesteps[0::indexSkipValue]
    timeToIndexArray=np.floor(1.0/dt*timestepsErr)
    timeToIndexArray = timeToIndexArray.astype(int)

    timestepsErrLog=np.logspace(2.0,np.log10(timesteps[-indexSkipValue]) , num=errorbarCount)
    timeToIndexArrayLog=np.floor(1.0/dt*timestepsErrLog)
    timeToIndexArrayLog=timeToIndexArrayLog.astype(int)

    for file in resultList:
        results = np.load(file)
        if results['maxEnergyError']>0.3:
            resultList.remove(file)
    
    stdMat = np.zeros((len(timeToIndexArray),len(resultList)))
    stdMatK = np.zeros((len(timeToIndexArray),len(resultList)))
    i=0
    for file in resultList:
            results = np.load(file)
            stdMat[:,i] = results['squaredQ'][timeToIndexArray] - results['aveQ'][timeToIndexArray]
            stdMatK[:,i] = results['K'][timeToIndexArray] 
    
            varQ += results['squaredQ'] - results['aveQ'] #not normalized yet
            varP += results['squaredP'] - results['aveP']
            K    = results['K']
            print(results['maxEnergyError'])
            i+=1
            print(i)
    
    std = np.std(stdMat, axis=1)
    stdK = np.std(stdMatK, axis=1)
#    print(np.sum(K))
    K*=1.0/np.sum(K)
    #K=moving_average(K,1000)    
    norm=1.0/(float(len(resultList)))
    varQ *= norm
    varP *= norm
    std  *= norm
    stdK  *= norm

    np.savez("../data/data", varQ=varQ, timesteps=timesteps, std=std, stdK=stdK, varP=varP, K=K, t1=t1, dt=dt, gamma=gamma)


else: 
    
    datafile=glob.glob('/users/stud/ledwon/Seafile/Aktuell/Masterarbeit/data/*.npz')
    data=np.load(datafile[0])
    varQ = data['varQ']
    timesteps=data['timesteps']
    std  = data['std']
    varP = data['varP']
    K    = data['K']
    stdK = data['stdK']
    t1=data['t1']
    dt=data['dt']
    gamma=data['gamma']
    
#linear 
    errorbarCount = 100
    indexSkipValue = int(1.0/dt * t1/float(errorbarCount))
    timestepsErr=timesteps[0::indexSkipValue]
    timeToIndexArray=np.floor(1.0/dt*timestepsErr)
    timeToIndexArray=timeToIndexArray.astype(int)

#log
    timestepsErrLog=np.logspace(2.0,np.log10(timesteps[-indexSkipValue]) , num=errorbarCount)
    timeToIndexArrayLog=np.floor(1.0/dt*timestepsErrLog)
    timeToIndexArrayLog=timeToIndexArrayLog.astype(int)




if gamma>1.0:
    diffType='super'
else:
    diffType='sub'

def memoryKernel(times):
    if diffType == 'super':
        result = np.cos((gamma)*np.arctan(times))*np.power(np.power(times,2)+1,-gamma/2.0)
        result *= -1.0/np.sum(result)
        return result

    if diffType == 'sub':
        result = np.power(times,-gamma)*np.power(scipy.special.gamma(gamma+1)*np.cos(0.5*np.pi*(gamma+1)),-1)
        result *= 1.0/np.sum(result)
        return result

#def theoDiff(times, gamma, fitindex):
#        return np.power(times,gamma) * varQ[fitindex]/np.power(times[fitindex],gamma)

def theoDiff(x,a,c):
    return a*np.power(x,gamma)+c

def linDiff(x,a,c):
    return a*x+c


startindex = int(math.floor((t1/dt)*0.20))
plotindex = int(math.floor((t1/dt)*0.05))
kernelplotindex = int(math.floor((t1/dt)*0.0))
endindex = int(math.floor(t1/dt)*0.7)
linindex = int(math.floor(t1/dt)*0.8)
popt, pcov = curve_fit(theoDiff, timesteps[startindex:endindex:2000],varQ[startindex:endindex:2000])
linpopt, linpcov = curve_fit(linDiff, timesteps[endindex::2000],varQ[endindex::2000])
print(popt)


var = plt.figure(1)
plt.xscale('log', nonposx="clip")
plt.yscale('log', nonposy="clip")
plt.plot(timesteps[plotindex::8000],varQ[plotindex::8000],label='Numerical results', color='#FC9169')
plt.errorbar(timestepsErrLog, varQ[timeToIndexArrayLog],yerr=std, fmt='none',capsize=1.0,ecolor='#FC9169')
plt.plot(timesteps[startindex:endindex:8000],theoDiff(timesteps[startindex:endindex:8000],popt[0],popt[1]), color='#0066FF',linestyle='--',label=r'$\propto t^{1.5}$')
plt.errorbar(timesteps[linindex::80000],linDiff(timesteps[linindex::80000],linpopt[0],linpopt[1]),linestyle=':',color='#009900',label=r'$\propto t$')
plt.xlabel('t')
plt.ylabel('Var(Q)')
plt.legend()
var.savefig("./img/varQlog.pdf",bbox_inches='tight')

var = plt.figure(2)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0), useMathText=True)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0), useMathText=True)
plt.plot(timesteps[::8000],varQ[::8000],label='Numerical results',color='#FC9169' )
plt.errorbar(timestepsErr, varQ[timeToIndexArray], yerr=std, fmt='none',capsize=1.0 ,ecolor='#FC9169')
plt.plot(timesteps[startindex:endindex:8000],theoDiff(timesteps[startindex:endindex:8000],popt[0],popt[1]),label=r'$\propto t^{1.5}$',color='#0066FF', linestyle='--')
plt.plot(timesteps[linindex::80000],linDiff(timesteps[linindex::80000],linpopt[0],linpopt[1]),label=r'$\propto t$', linestyle=':',color='#009900')
plt.xlabel('t')
plt.ylabel('Var(Q)')
plt.legend()
var.savefig("./img/varQ.pdf",bbox_inches='tight')


kerneltimes=np.logspace(timesteps[kernelplotindex],timesteps[-1],50000)
Kernel = plt.figure(3)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0), useMathText=True)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0), useMathText=True)
plt.xscale('log', nonposx="clip")
plt.yscale('log', nonposy="clip")
#plt.plot(timestepsErr,np.abs(K[timeToIndexArray]),label='Bath memory kernel',color='#FC9169' )
plt.plot(timestepsErr,K[timeToIndexArray],label='Bath memory kernel',color='#FC9169' )
plt.errorbar(timestepsErr, np.abs(K[timeToIndexArray]),yerr=stdK, fmt='none',capsize=1.0,ecolor='#FC9169',elinewidth='0.7')
plt.plot(kerneltimes,np.abs(memoryKernel(kerneltimes)),label='Theoretical memory kernel', linestyle=':')
plt.xlabel('t')
plt.ylabel('Memory Kernel')
plt.legend()
Kernel.savefig("./img/K.pdf",bbox_inches='tight')
print(K[0])

#plt.show()
