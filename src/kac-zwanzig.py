import numpy as np
import matplotlib.pyplot as plt
import math
from cmath import *
from heatbath import *
from bathensemble import *
import scipy


########################################################################################################################
#set parameters
n=5#ensemble of baths that we average over
N=8000 #number of bath oscillators
beta=1.0 #1\kB*T
Q0=0.0 #starting pos/impulse of distinguished particle
P0=0.0
oscMass=1.0 #1.0 #mass of heaviest bath oscillator
M=0.01# mass of the distinguished particle
#masses=m*np.ones(N)
t0=0.1
t1=3000.0
dt=0.0005#3.0/float(N)#(t1-t0)/100.0


gamma=1.3

if gamma>1.0:
    diffType='super'
else:
    diffType='sub'


startstr='Starting Simulation, averaging over ' + str(n) + ' heatbaths with '+ str(N) + ' oscillators'
print(startstr)


ensemble1 = bathensemble(n,N,beta,Q0,P0,oscMass,M,t0,t1,dt,gamma,diffType)
ensemble1.averageEnsemble()

#gamma=1.4
#diffType='super' #'super' for superdiffusion, 'sub' for subdiffusion
#ensemble2 = bathensemble(n,N,beta,Q0,P0,oscMass,M,t0,t1,dt,gamma,diffType)
#ensemble2.averageEnsemble()


#np.save('varQarray',ensemble1.varQ)

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
        return np.power(times,gamma) * ensemble1.varQ[math.floor((t1/dt)*0.3)]/np.power(ensemble1.timesteps[math.floor((t1/dt)*0.3)],gamma)


#startindex = math.floor((t1/dt)*0.2)
#const = ensemble1.varQ[math.floor((t1/dt)*0.3)]/np.power(ensemble1.timesteps[math.floor((t1/dt)*0.3)],gamma)
#const_control = ensemble1.varQ[math.floor((t1/dt)*0.3)]/np.power(ensemble1.timesteps[math.floor((t1/dt)*0.3)],1.5)
#plt.figure(1)
#plt.loglog(ensemble1.timesteps,ensemble1.varQ,ensemble1.timesteps[math.floor((t1/dt)*0.1):math.floor(t1/dt)],theoDiff(ensemble1.timesteps,gamma)[math.floor((t1/dt)*0.1):math.floor(t1/dt)],ensemble1.timesteps[math.floor((t1/dt)*0.1):math.floor(t1/dt)],theoDiff(ensemble1.timesteps,1.0)[math.floor((t1/dt)*0.1):math.floor(t1/dt)])
##plt.loglog(ensemble1.timesteps,ensemble1.varQ,ensemble1.timesteps[math.floor((t1/dt)*0.1):math.floor(t1/dt)],theoDiff(ensemble1.timesteps,const,gamma)[math.floor((t1/dt)*0.1):math.floor(t1/dt)],ensemble1.timesteps[math.floor((t1/dt)*0.1):math.floor(t1/dt)])#,theoDiff(ensemble1.timesteps,const,1.4)[math.floor((t1/dt)*0.1):math.floor(t1/dt)])
##plt.hist(ensemble.singleBath.q[0,:],1000)
#plt.xlabel('t')
#plt.ylabel('Var(Q)')
#
#
#plt.figure(2)
#plt.plot(ensemble1.timesteps,ensemble1.K,ensemble1.timesteps,memoryKernel(ensemble1.timesteps))
##plt.hist(ensemble.singleBath.q[0,:],1000)
#plt.xlabel('t')
#plt.ylabel('Memory Kernel')
#
#plt.figure(3)
#plt.plot(ensemble1.timesteps,ensemble1.singleBath.Q)
#
#plt.figure(4)
#plt.plot(ensemble1.timesteps,ensemble1.singleBath.energy)
#
#plt.figure(5)
#plt.plot(ensemble1.timesteps,ensemble1.singleBath.momentum)
#
#plt.show()












#to do/ideas:
#replicate Subdiffusion from kupfermann paper
#change omega distribution from uniform to something else
#change values of c & a
#real f function
