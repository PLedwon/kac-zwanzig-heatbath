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
N=16000 #number of bath oscillators
beta=1.0 #1\kB*T
Q0=0.0 #starting pos/impulse of distinguished particle
P0=0.0
oscMass=1.0 #1.0 #mass of heaviest bath oscillator
M=0.01# mass of the distinguished particle
#masses=m*np.ones(N)
t0=0.1
t1=3000.0
dt=0.0005#3.0/float(N)#(t1-t0)/100.0


gamma=1.5

if gamma>1.0:
    diffType='super'
else:
    diffType='sub'


startstr='Starting Simulation, averaging over ' + str(n) + ' heatbaths with '+ str(N) + ' oscillators'
print(startstr)


ensemble1 = bathensemble(n,N,beta,Q0,P0,oscMass,M,t0,t1,dt,gamma,diffType)
ensemble1.averageEnsemble()

