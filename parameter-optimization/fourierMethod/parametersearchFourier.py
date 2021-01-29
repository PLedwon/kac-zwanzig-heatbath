import numpy as np
import matplotlib.pyplot as plt
import math
from cmath import *
import scipy.special
from scipy.fftpack import fft,ifft


########################################################################################################################
#set parameters
N=20000 #number of bath oscillators
beta=1.0 #1\kB*T
Q0=0.0 #starting pos/impulse of distinguished particle
P0=0.0
oscMass=1.0 #1.0 #mass of heaviest bath oscillator
M=0.01# mass of the distinguished particle
#masses=m*np.ones(N)
t0=0.0
t1=5000.0
dt=2.0#01.0/float(N)#(t1-t0)/100.0
dt=(t1-t0)/float(N)
Omega=1.0
gridsize = 51 #should be (M*10)-1 for nice values
timesteps=np.arange(t0,t1,dt)
#timesteps=np.logspace(0.0,np.log10(t1),2000)#np.arange(0,t1,1.0)
lowerNRange = np.linspace(-1.5,-1.0,gridsize)
upperNRange = np.linspace(1.2,1.5,gridsize)
#lowerNRange =np.arange(-0.882,-0.880,0.0005)
#upperNRange =np.arange(1.171,1.173,0.0005)
cutoff = 10000
kernelDiff = cutoff*np.ones((len(lowerNRange),len(upperNRange)))

gamma=1.2


def memoryKernel(timesteps):
        result = np.cos((gamma)*np.arctan(timesteps*Omega))*np.power(np.power(timesteps/Omega,2)+1,-gamma/2.0)
        result *= 1.0/np.sum(result)
        return result

massKifft=ifft(memoryKernel(np.arange(t0,t1,(t1-t0)/float(N))))
realK = memoryKernel(timesteps)
def computeMassesFourierFourier(omega):
       return oscMass * np.power(omega,-2.0) * massKifft

def computeMasses(omega)

def setFrequencyRange(a,b):
            omega_min=N**a
            omega_max=omega_min*N**b
            return omega_min, omega_max

def computeKernel(timesteps,k,omega):
        K = np.zeros(timesteps.size)
        for i in range(0,timesteps.size):
            K[i]= np.dot(k,np.cos(omega*timesteps[i]))
        K *= 1.0/np.sum(K)
        return K


invK = np.reciprocal(realK)

counter=0

for i in range(0,len(lowerNRange)-1):
    for j in range(0,len(upperNRange)-1):

            a=lowerNRange[i]
            b=upperNRange[j]
            omega_min, omega_max = setFrequencyRange(a,b)
            omega = np.linspace(omega_min,omega_max,num=N)
            #omega =np.random.uniform(omega_min,omega_max,N) # np.linspace(omega_min,omega_max,num=N)
            masses = computeMassesFourier(omega)
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
masses = computeMassesFourier(omega)
timesteps=np.arange(t0,t1,1.0)
k=np.multiply(masses,np.power(omega,2)) # compute spring constants
K = computeKernel(timesteps,k,omega)
realK = memoryKernel(timesteps)

#np.savez('parameters', ind=ind, kernelDiff=kernelDiff, lowerNRange=lowerNRange, upperNRange=upperNRange)

kern = plt.figure(1)
plt.plot(timesteps,K,timesteps,realK)
plt.xlabel('t')
plt.ylabel('Memory Kernel')
kern.savefig("MemoryKernel.pdf")

timestepsLog=np.logspace(0.0,np.log10(timesteps[-1]),2000)#np.arange(0,t1,1.0)



K = computeKernel(timestepsLog,k,omega)
realK = memoryKernel(timestepsLog)
kernlog = plt.figure(2)
plt.plot(timestepsLog,np.abs(K),timestepsLog,np.abs(realK))
plt.xlabel('t')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Memory Kernel')
kernlog.savefig("MemoryKernelLog.pdf")

mat=plt.figure(8)
plt.imshow(kernelDiff)
plt.colorbar
#plt.show()

diff = plt.figure(7)
plt.plot(timestepsLog,np.abs(K-realK)/np.sum(realK),label=r'$\Delta K/K, \gamma=1.2$')
plt.xlabel('t')
plt.ylabel('rel. Error')
plt.xscale('log')
plt.yscale('log')
diff.savefig("relErr.pdf")
