import numpy as np
import glob
import matplotlib.pyplot as plt
import math
import scipy

resultList = glob.glob('*.np[yz]')

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

for file in resultList:
    results = np.load(file)
    varQ += results['squaredQ'] - results['aveQ'] #have to normalize
    varP += results['squaredP'] - results['aveP']
    K    += results['K']

K*=1.0/np.sum(K)
#parallelizedn = 6
norm=1.0/float(len(results))
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

def theoDiff(times,gamma):
        return np.power(times,gamma) * varQ[math.floor((t1/dt)*0.3)]/np.power(timesteps[math.floor((t1/dt)*0.3)],gamma)



startindex = math.floor((t1/dt)*0.2)
const = varQ[math.floor((t1/dt)*0.3)]/np.power(timesteps[math.floor((t1/dt)*0.3)],gamma)
#const_control = varQ[math.floor((t1/dt)*0.3)]/np.power(ensemble1.timesteps[math.floor((t1/dt)*0.3)],1.5)
var = plt.figure(1)
plt.plot(timesteps,varQ,timesteps[math.floor((t1/dt)*0.1):math.floor(t1/dt)],theoDiff(timesteps,gamma)[math.floor((t1/dt)*0.1):math.floor(t1/dt)])
#plt.loglog(ensemble1.timesteps,ensemble1.varQ,ensemble1.timesteps[math.floor((t1/dt)*0.1):math.floor(t1/dt)],theoDiff(ensemble1.timesteps,const,gamma)[math.floor((t1/dt)*0.1):math.floor(t1/dt)],ensemble1.timesteps[math.floor((t1/dt)*0.1):math.floor(t1/dt)])#,theoDiff(ensemble1.timesteps,const,1.4)[math.floor((t1/dt)*0.1):math.floor(t1/dt)])
#plt.hist(ensemble.singleBath.q[0,:],1000)
plt.xlabel('t')
plt.ylabel('Var(Q)')
var.savefig("varQ.pdf",bbox_inches='tight')


K = plt.figure(2)
plt.plot(timesteps,memoryKernel(timesteps))
plt.xlabel('t')
plt.ylabel('Memory Kernel')
K.savefig("K.pdf",bbox_inches='tight')

#plt.figure(3)
#plt.plot(timesteps,ensemble1.singleBath.Q)
#
#plt.figure(4)
#plt.plot(ensemble1.timesteps,ensemble1.singleBath.energy)
#
#plt.figure(5)
#plt.plot(ensemble1.timesteps,ensemble1.singleBath.momentum)

#plt.show()
