import numpy as np
import matplotlib.pyplot as plt
import math
from cmath import *
import scipy.special


########################################################################################################################
#set parameters
N=16000 #number of bath oscillators
beta=1.0 #1\kB*T
Q0=0.0 #starting pos/impulse of distinguished particle
P0=0.0
oscMass=1.0 #1.0 #mass of heaviest bath oscillator
M=0.01# mass of the distinguished particle
#masses=m*np.ones(N)
t0=0.1
t1=3000.0
dt=0.001#01.0/float(N)#(t1-t0)/100.0
#dt=5.5
gridsize = 1
timesteps=np.arange(0.0,t1,dt)
lowerNRange = np.linspace(-1.1,-0.9,gridsize)
upperNRange = np.linspace(1.0,1.3,gridsize)
#lowerNRange =np.arange(-1.1,-0.9,0.1)
#upperNRange =np.arange(0.8,1.3,0.1)
cutoff = 10000
kernelDiff = cutoff*np.ones((len(lowerNRange),len(upperNRange)))

gamma=1.5

if gamma>1.0:
    diffType='super'
else:
    diffType='sub'

def computeMasses(omega):
  if diffType == 'sub':
       return oscMass * np.power(omega/np.amin(omega),gamma-3.0)
  if diffType == 'super':
       Omega = 1.0
       return oscMass * np.power(omega/np.amin(omega),gamma-3.0)*np.exp(-omega/Omega)

def setFrequencyRange(a,b):
        if diffType == 'sub':
            omega_min=N**a
            omega_max=omega_min*N**b

        if diffType == 'super':
            omega_min=N**a
            omega_max=omega_min*N**b

        return omega_min, omega_max

def computeKernel(timesteps,k,omega):
        K = np.zeros(timesteps.size)
        for i in range(0,timesteps.size):
            K[i]= np.sum(np.multiply(k,np.cos(omega*timesteps[i])))
        K *= 1.0/np.sum(K)
        return K

def memoryKernel(timesteps):
    if diffType == 'super':
        result = np.cos((gamma)*np.arctan(timesteps))*np.power(np.power(timesteps,2)+1,-gamma/2.0)
        result *= 1.0/np.sum(result)
        return result

    if diffType == 'sub':
        result = np.power(timesteps,-gamma)*np.power(scipy.special.gamma(gamma+1)*np.cos(0.5*np.pi*(gamma+1)),-1)
        result *= 1.0/np.sum(result)
        return result

realK = memoryKernel(timesteps)
invK = np.reciprocal(realK)

counter=0

for i in range(0,len(lowerNRange)-1):
    for j in range(0,len(upperNRange)-1):

            a=lowerNRange[i]
            b=upperNRange[j]
            omega_min, omega_max = setFrequencyRange(a,b)
            omega = np.linspace(omega_min,omega_max,num=N)
            #omega =np.random.uniform(omega_min,omega_max,N) # np.linspace(omega_min,omega_max,num=N)
            masses = computeMasses(omega)
            k=np.multiply(masses,np.power(omega,2)) # compute spring constants
            K = computeKernel(timesteps,k,omega)
            kernelDiff[i,j] = dt*np.sum(np.abs(K-realK))/np.sum(realK)/t1
            if kernelDiff[i,j] > cutoff:
                    kernelDiff[i,j] =cutoff
            counter+=1
            progressstr= str(counter) + '/'+ str(np.power(gridsize,2)) + ' done'
            print(progressstr)
ind = np.unravel_index(np.argmin(kernelDiff, axis=None), kernelDiff.shape)
print(kernelDiff[ind])
print(lowerNRange[ind[0]],upperNRange[ind[1]])


omega_min, omega_max = setFrequencyRange(lowerNRange[ind[0]],upperNRange[ind[1]])
omega =np.linspace(omega_min,omega_max,num=N)
masses = computeMasses(omega)
timesteps=np.arange(t0,t1,0.1)
k=np.multiply(masses,np.power(omega,2)) # compute spring constants
K = computeKernel(timesteps,k,omega)
realK = memoryKernel(timesteps)

kern = plt.figure(1)
plt.plot(timesteps,K,timesteps,realK)
plt.xlabel('t')
plt.ylabel('Memory Kernel')
kern.savefig("MemoryKernel.pdf")

timesteps=np.arange(0,2*t1,0.1)
K = computeKernel(timesteps,k,omega)
realK = memoryKernel(timesteps)
kernlog = plt.figure(2)
plt.plot(timesteps,np.abs(K),timesteps,np.abs(realK))
plt.xlabel('t')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Memory Kernel')
kernlog.savefig("MemoryKernelLog.pdf")
#plt.show()
