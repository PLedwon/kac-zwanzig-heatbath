import numpy as np
import matplotlib.pyplot as plt
import math
from cmath import *
from heatbath import *
from bathensemble import *
import scipy


########################################################################################################################
#set parameters
n=1#ensemble of baths that we average over
N=2 #number of bath oscillators
beta=1.0 #1\kB*T
M=np.power(10.0,-2)# mass of the distinguished particle
Q0=0 #starting pos/impulse of distinguished particle
P0=np.power(beta,-0.5)*np.power(M,0.5)*np.random.standard_normal(1)
oscMass=np.power(10.0,2) #1.0 #mass of heaviest bath oscillator
#masses=m*np.ones(N)
t0=0.0
t1=np.power(10.0,0.0)
dt=np.power(10.0,-5.0)#0.0004#3.0/float(N)#(t1-t0)/100.0
Omega=1.0
omega_min=N**(-0.8323)
omega_max=omega_min*N**1.05764

gamma=1.2

if gamma>1.0:
    diffType='super'
else:
    diffType='sub'


startstr='Starting Simulation, averaging over ' + str(n) + ' heatbaths with '+ str(N) + ' oscillators'
print(startstr)


ensemble1 = bathensemble(n,N,beta,Q0,P0,oscMass,M,t0,t1,dt,gamma,diffType,Omega,omega_min,omega_max)
ensemble1.averageEnsemble()
print(ensemble1.singleBath.Q)
