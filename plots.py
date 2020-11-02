#!/usr/bin/python3
import numpy as np
import glob
import matplotlib.pyplot as plt
import math
import scipy
from scipy.optimize import curve_fit

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
timestepsErr=timesteps[::indexSkipValue]
timeToIndexArray=np.floor(1.0/dt*timestepsErr)
timeToIndexArray = timeToIndexArray.astype(int)
for file in resultList:
    results = np.load(file)
    if results['maxEnergyError']>0.3:
        resultList.remove(file)

stdMat = np.zeros((len(timeToIndexArray),len(resultList)))
i=0
for file in resultList:
        results = np.load(file)
        stdMat[:,i] = results['squaredQ'][timeToIndexArray] - results['aveQ'][timeToIndexArray]

        varQ += results['squaredQ'] - results['aveQ'] #not normalized yet
        varP += results['squaredP'] - results['aveP']
        K    += results['K']
        print(results['maxEnergyError'])
        i+=1
        print(i)

std = np.std(stdMat, axis=1)
K*=1.0/np.sum(K)
norm=1.0/(float(len(resultList)))
print(len(resultList))
varQ *= norm
varP *= norm
std  *= norm


if gamma>1.0:
    diffType='super'
else:
    diffType='sub'

def memoryKernel(times):
    if diffType == 'super':
        result = np.cos((gamma)*np.arctan(times))*np.power(np.power(times,2)+1,-gamma/2.0)
        result *= 1.0/np.sum(result)
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


startindex = int(math.floor((t1/dt)*0.00))
endindex = int(math.floor(t1/dt)*0.7)
linindex = int(math.floor(t1/dt)*0.8)
#fitindex = int(math.floor((t1/dt)*0.2))
popt, pcov = curve_fit(theoDiff, timesteps[startindex:endindex:2000],varQ[startindex:endindex:2000])
linpopt, linpcov = curve_fit(linDiff, timesteps[endindex::2000],varQ[endindex::2000])
print(popt)

var = plt.figure(1)
plt.loglog(timesteps[startindex::8000],varQ[startindex::8000],label='Numerical results')
plt.errorbar(timestepsErr, varQ[timeToIndexArray], std, fmt='none', elinewidth='0.2')
plt.loglog(timesteps[startindex:endindex:8000],theoDiff(timesteps[startindex:endindex:8000],popt[0],popt[1]), linestyle='--',label=r'$\propto t^{1.5}$')
plt.loglog(timesteps[linindex::80000],linDiff(timesteps[linindex::80000],linpopt[0],linpopt[1]),linestyle=':',color='r',label=r'$\propto t$')
plt.xlabel('t')
plt.ylabel('Var(Q)')
plt.legend()
var.savefig("./img/varQlog.pdf",bbox_inches='tight')

var = plt.figure(2)
plt.plot(timesteps[::8000],varQ[::8000],label='Numerical results',color='#FC9169' )
plt.errorbar(timestepsErr, varQ[timeToIndexArray], std, fmt='none', ecolor='#FC9169',elinewidth='0.2')
#plt.plot(timesteps[startindex:endindex:8000],theoDiff(timesteps,gamma,fitindex)[startindex:endindex:8000],label=r'$\propto t^{1.5}$')
plt.plot(timesteps[startindex:endindex:8000],theoDiff(timesteps[startindex:endindex:8000],popt[0],popt[1]),label=r'$\propto t^{1.5}$', linestyle='--')
plt.plot(timesteps[linindex::80000],linDiff(timesteps[linindex::80000],linpopt[0],linpopt[1]),label=r'$\propto t$', linestyle=':',color='r')
plt.xlabel('t')
plt.ylabel('Var(Q)')
plt.legend()
var.savefig("./img/varQ.pdf",bbox_inches='tight')

Kernel = plt.figure(3)
plt.plot(timesteps[::1600],K[::1600],label='Bath memory kernel')
plt.plot(timesteps[::16000],memoryKernel(timesteps)[::16000],label='Theoretical memory kernel', linestyle=':')
plt.xlabel('t')
plt.ylabel('Memory Kernel')
plt.legend()
Kernel.savefig("./img/K.pdf",bbox_inches='tight')


#plt.show()
