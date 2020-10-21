#!/usr/bin/python3
import numpy as np
import glob
import matplotlib.pyplot as plt
import math
import scipy

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
errorFiles=0

for file in resultList:
    results = np.load(file)
    if results['maxEnergyError']>0.1:
        errorFiles+=1
    
    else:
        varQ += results['squaredQ'] - results['aveQ'] #have to normalize
        varP += results['squaredP'] - results['aveP']
        K    += results['K']
        print(results['maxEnergyError'])

K*=1.0/np.sum(K)
norm=1.0/(float(len(results))-errorFiles)
varQ *= norm
varP *= norm


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

def theoDiff(times, gamma, fitindex):
        return np.power(times,gamma) * varQ[fitindex]/np.power(times[fitindex],gamma)



startindex = int(math.floor((t1/dt)*0.10))
endindex = int(math.floor(t1/dt)*0.85)
fitindex = int(math.floor((t1/dt)*0.50))

var = plt.figure(1)
plt.loglog(timesteps[startindex:endindex:16000],varQ[startindex:endindex:16000],':',timesteps[startindex:endindex:16000],theoDiff(timesteps,gamma,fitindex)[startindex:endindex:16000])
plt.xlabel('t')
plt.ylabel('Var(Q)')
var.savefig("./img/varQlog.pdf",bbox_inches='tight')

var = plt.figure(2)
plt.plot(timesteps[::16000],varQ[::16000],':',timesteps[startindex:endindex:16000],theoDiff(timesteps,gamma,fitindex)[startindex:endindex:16000])
plt.xlabel('t')
plt.ylabel('Var(Q)')
var.savefig("./img/varQ.pdf",bbox_inches='tight')

Kernel = plt.figure(3)
plt.plot(timesteps[::16000],memoryKernel(timesteps)[::16000],timesteps[::1600],K[::1600],'.')
plt.xlabel('t')
plt.ylabel('Memory Kernel')
Kernel.savefig("./img/K.pdf",bbox_inches='tight')


#plt.show()
