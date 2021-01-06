import numpy as np
import matplotlib.pyplot as plt
from cmath import *
from heatbath import *

class bathensemble():
    def __init__(self,n,N,beta,Q0,P0,oscMass,M,t0,t1,dt,gamma,diffType,Omega,omega_min,omega_max):

        self.n = n # number of baths in the ensemble
        self.N=N #number of bath oscillators
        self.beta=beta #1\kB*T
        self.Q0=Q0 #starting pos/impulse of distinguished particle
        self.P0=P0
        self.oscMass=oscMass #mass of the heaviest bath oscillator
        self.M=M
        self.t0=t0
        self.t1=t1
        self.dt=dt
        self.gamma=gamma
        self.timesteps=np.arange(self.t0,self.t1,self.dt)
        self.diffType = diffType
        self.Omega = Omega
        self.omega_min=omega_min
        self.omega_max=omega_max
        self.K_N=np.zeros(self.timesteps.size) # memory kernel value at every timestep
        self.K = np.zeros(self.timesteps.size)
        self.singleBathK_N = np.zeros(self.timesteps.size)
        


        self.aveQ=np.zeros(self.timesteps.size)
        self.aveP=np.zeros(self.timesteps.size)
        self.squaredQ=np.zeros(self.timesteps.size)
        self.squaredP=np.zeros(self.timesteps.size)
        self.varQ=np.zeros(self.timesteps.size)
        self.varP=np.zeros(self.timesteps.size)


    def computeMasses(self,omega):
        if self.diffType == 'sub':
            return self.oscMass * np.power(omega/np.amin(omega),self.gamma-3.0)#np.power(omega,1.0-self.alpha)*self.N**(self.a-1)

        if self.diffType == 'super':
            return self.oscMass * np.power(omega/np.amin(omega),self.gamma-3.0)*np.exp(-omega/self.Omega)#np.power(omega,self.alpha)*np.exp(-omega*1.0/Omega)*self.N**(self.a-1)


    def simulateSingleBath(self):

        omega =np.linspace(self.omega_min,self.omega_max,num=self.N)
        self.omega=omega
        #omega =np.random.uniform(omega_min,omega_max,self.N) 

        masses = self.computeMasses(omega)
        k=np.multiply(masses,np.power(omega,2)) # compute spring constants
        self.k=k
        #fig=plt.figure(3)
        #plt.hist(k,200)

        #generate starting positions/impulses for oscillators
        q0 = self.Q0 + np.power(self.beta,-0.5)*np.multiply(np.power(k,-0.5),np.random.standard_normal(self.N))
        p0 = np.power(self.beta,-0.5)*np.multiply(np.power(masses,0.5),np.random.standard_normal(self.N))
        p0 -= np.average(p0)

        # cast initial values into one array
        y0=np.hstack([self.Q0,self.P0,q0,p0])
        self.singleBath = heatbath(self.N,y0,k,masses,self.M,self.t0,self.t1,self.dt)

       # self.singleBathK_N = np.zeros(self.timesteps.size)
       # for i in range(0,self.timesteps.size):
       #     self.singleBathK_N[i]= np.dot(k,np.cos(omega*self.timesteps[i]))
       # self.singleBathK_N *= 1.0/np.sum(self.singleBathK_N)


    def averageEnsemble(self):

        self.maxEnergyError = 0
        self.maxMomentumError = 0
        self.avgEnergyError = 0
        for i in range(0,self.n):
            print(i)
            self.simulateSingleBath()
            self.aveQ += 1.0/float(self.n) *self.singleBath.Q
            self.aveP += 1.0/float(self.n) *self.singleBath.P
            self.squaredQ += 1.0/float(self.n) *np.power(self.singleBath.Q,2)
            self.squaredP += 1.0/float(self.n) *np.power(self.singleBath.P,2)

            #self.K = self.singleBathK_N
            #self.K += 1.0/float(self.n) * self.singleBathK_N
            self.avgEnergyError +=1.0/float(self.n) *self.singleBath.avgEnergyError

            if self.maxEnergyError<self.singleBath.maxEnergyError:
                self.maxEnergyError=self.singleBath.maxEnergyError

            if self.maxMomentumError<self.singleBath.maxMomentumError:
                self.maxMomentumError=self.singleBath.maxMomentumError


        #self.varQ = self.squaredQ - np.power(self.aveQ,2)
        #self.varP = self.squaredP - np.power(self.aveP,2)
        #print('max energy error = ' + str(math.ceil(self.maxEnergyError*10000.0)/100.0) + '%' )
        #print('avg energy error = ' + str(math.ceil(self.avgEnergyError*10000.0)/100.0) + '%' )
        #print('max momentum error: ' + str(self.maxMomentumError))

        print('max energy error =' )
        print(self.maxEnergyError)
        print('avg energy error = ' )
        print(self.avgEnergyError)
        print(self.maxMomentumError)

        name = str(np.floor(np.random.uniform(0,999999,1)))
        np.savez(name, squaredQ=self.squaredQ, squaredP=self.squaredP, aveQ=self.aveQ, aveP=self.aveP, maxEnergyError=self.maxEnergyError, maxMomentumError=self.maxMomentumError, dt=self.dt, t1=self.t1,timesteps=self.timesteps, gamma=self.gamma, avgEnergyError=self.avgEnergyError,Omega=self.Omega, omega=self.omega,k=self.k)
