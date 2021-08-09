import numpy as np
import math
from cmath import *
from scipy.integrate import odeint


#No Potential yet

class heatbath():
    def __init__(self,N,y0,k,masses,M,t0,t1,dt):
        self.N=N                                    # number of oscillators the bath consists of
        self.k=k                    # vector of spring constant of the oscillators
        self.masses=masses                    # vector of masses of the oscillators, must not be int, reciprocal doenst work with ints
        self.invm=np.reciprocal(self.masses)  #invert masses for maketimesteps-function, requires m to be floatarray
        self.invM=1.0/M
        self.t0=t0
        self.t1=t1
        self.y0=y0                          # y0=(Q0,P0,q0,p0)
        self.dt=dt
        self.timesteps=np.arange(self.t0,self.t1,self.dt)
        self.initialEnergy = 0.5*(self.invM*np.power(self.y0[0],2)+np.inner(self.invm,np.power(self.y0[self.N+2:2*self.N+2],2)) +  np.inner(self.k,np.power((self.y0[2:self.N+2]-self.y0[0]),2)))     # Initial energy of the system
        self.energy= np.zeros(len(self.timesteps))
        self.momentum = np.zeros(len(self.timesteps))
        self.dydt=np.zeros(2*self.N+2)

        def runHeatbath(self):
            self.solve_ivp()
            self.checkEnergy()
            self.checkMomentum()

        runHeatbath(self)

##########################################################
#    def hamiltonEOM(self,y,t):
#
#        self.dydt[0]                   = y[1]*self.invM
#        self.dydt[1]                   = np.inner(self.k,(y[2:self.N+2]-y[0])) #+0.5- y[0]# - np.power(y[0],3) # + potential can be added
#        self.dydt[2:self.N+2]          = np.multiply(self.invm,y[self.N+2:2*self.N+2])
#        self.dydt[self.N+2:2*self.N+2] = -np.multiply(self.k,(y[2:self.N+2]-y[0]))
#        return self.dydt
#
#    def solve_ivp(self):
#        self.sol = odeint(self.hamiltonEOM,self.y0,self.timesteps)
#        self.Q=self.sol[:,0]
#        self.P=self.sol[:,1]
#        self.q=self.sol[:,2:self.N+2]
#        self.p=self.sol[:,self.N+2:2*self.N+2]

###########################################################

#symplectic euler version
    def solve_ivp(self):

        self.q=np.zeros((self.N,1))
        self.p=np.zeros((self.N,1))
        self.Q=np.zeros(len(self.timesteps))
        self.P=np.zeros(len(self.timesteps))

        self.q=self.y0[1:self.N+1]
        self.p=self.y0[self.N+2:2*self.N+2]
        self.Q[0]=self.y0[0]
        self.P[0]=self.y0[self.N+1]
        print(self.P[0])

        self.energy[0]=self.initialEnergy
        self.momentum[0]=np.sum(self.p)+self.P[0]
        #self.energyE
    #semi-implicit Euler-scheme, might be replaced by an other integrator
    #impulses get updated first
        for i in range(0,len(self.timesteps)-1): #Hamilton eom in the euler scheme
            self.p = self.p - np.multiply(self.k,self.q-self.Q[i]) *self.dt
            self.P[i+1]   = self.P[i]   + np.inner(self.k,self.q-self.Q[i])    *self.dt
            self.Q[i+1]   = self.Q[i]   + self.P[i+1]                               *self.dt
            self.q = self.q + np.multiply(self.p,self.invm)      *self.dt

            #energies and momentum errors
            self.energy[i+1] = 0.5*(self.invM*np.power(self.P[i],2)+ np.inner(self.invm,np.power(self.p,2)) +  np.inner(self.k,np.power((self.q-self.Q[i]),2)))
            self.momentum[i+1] = np.sum(self.p) + self.P[i]

##############################################

    def checkEnergy(self): # error of the energy for every timestep
        self.energyError = np.reciprocal(self.initialEnergy)*np.abs(self.energy-self.initialEnergy)
        self.maxEnergyError=np.max(self.energyError)
        self.avgEnergyError=np.average(self.energyError)

    def checkMomentum(self):
        self.momentumError=np.abs(self.momentum-self.momentum[0])#/(self.P[0]-np.sum(self.p[0,:]))
        self.maxMomentumError=np.max(self.momentumError)
