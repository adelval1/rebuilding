import numpy as np
from numpy import random
import matplotlib.pyplot as plt

class metropolis:

    """

        This class creates a Metropolis-Hastings sampler (MCMC chain) with classical adaptation of the covariance matrix through a burn-in stage.

        INPUTS: Initial covariance matrix, the posterior (log-posterior, proportional), and number of samples to adapt the covariance matrix (burning stage)

        OUTPUTS: seed function to set some initial parameters of the chain, move one step by an accept-reject procedure, burning stage whose samples are used to adapt the covariance matrix. Other utilities are included to ease the use and bug finding.

    """

    def __init__(self, covinit, fun_in, nburn):
        self.cov = covinit
        self.fun = fun_in
        self.Tac = [0.]*nburn
        self.d = len(covinit)
        self.CovDev = np.linalg.cholesky(self.cov*(2.38*2.38)/self.d)
        # np.random.seed(2) # Comment or uncomment depending if you are trying to compare results and don't want stochastic noise.

    def seed(self,xini):
        self.xcur = xini
        self.fcur = self.fun(xini)
        self.mea = self.xcur
        self.acc = np.array([[0.]*len(self.xcur)]*len(self.xcur))
        for i in range(len(self.xcur)):
            for j in range(len(self.xcur)):
                self.acc[i][j] += (self.xcur[i]*self.xcur[j])
        self.count = 1

    def DoStep(self,nstep):
        for i in range(nstep):
            self.doOneStep()
        return self.xcur

    def doOneStep(self):
        xTemp = np.add(self.xcur, self.PropStep())
        fTemp = self.fun(xTemp)
        r0=np.random.random()
        if fTemp-self.fcur>np.log(r0):
            self.xcur = xTemp
            self.fcur = fTemp
            self.Tac[self.count%len(self.Tac)] = 1
        else:
            self.Tac[self.count%len(self.Tac)] = 0

    def Burn(self):
        print("Burning "+str(len(self.Tac))+" samples...")
        for i in range(len(self.Tac)):
            self.BurnOneStep()
            if i>300 and i%100==0: ## If iteration number is larger than 300 and a multiplo of 100
                self.AdaptCov()
        print("Accept rate " + str(np.sum(self.Tac)/len(self.Tac)) +"\n")

    def BurnOneStep(self):
        xTemp = np.add(self.xcur,self.PropStep())
        fTemp = self.fun(xTemp)
        r0 = np.random.random()
        if fTemp - self.fcur > np.log(r0):
            self.xcur = xTemp
            self.fcur = fTemp
            self.Tac[self.count%len(self.Tac)] = 1
        else:
            self.Tac[self.count%len(self.Tac)] = 0

        self.mea = np.add(self.mea,self.xcur)
        for i in range(len(self.xcur)):
            for j in range(len(self.xcur)):
                self.acc[i][j] += (self.xcur[i]*self.xcur[j])
        self.count += 1

    def GetCur(self):
        return self.xcur

    def SetCovProp(self,Cov_in):
        self.cov = Cov_in

    def GetCovProp(self):
        return self.cov

    def DecCov(self):
        self.CovDev = np.linalg.cholesky(self.cov*(2.38*2.38)/self.d)

    def PropStep(self):
        dX = [0.]*self.d
        for i in range(self.d):
            dX[i] = np.random.normal()

        for i in range(len(dX)):
            dX[i] = np.sum([self.CovDev[i,j]*dX[j] for j in range(len(dX))])

        # dX = np.matmul(self.CovDev,dX)
        for i in range(self.d):
            dX[i] += np.random.normal()*1.e-6
        return dX

    def AdaptCov(self):
        c1 = self.count
        c2 = c1*c1
        self.cov = self.acc/(self.count-1)
        for i in range(len(self.xcur)):
            for j in range(len(self.xcur)):
                self.cov[i][j] -= (self.mea[i]*self.mea[j])/c2
        self.cov += 1.e-6 * np.identity(self.d)
        self.DecCov()