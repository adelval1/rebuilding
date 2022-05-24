import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize
import sampler as mcmc
import subprocess
from subprocess import call
import os

def BL(H,Tw,Pd,Ps,gamma,start_from_previous):

    if start_from_previous == True:
    	os.system("scp ./neboulaoutput/stagsol.dat ./neboulainput")
    	os.system("mv ./neboulainput/stagsol.dat ./neboulainput/bc_ini.in")

    #---------------------------------
    # Write gamma for BL code on wall conditions

    with open('./neboulainput/bc_wall.in', 'r') as f:
        lines = f.readlines()

    # Replace the lines
    lines[3] = ('{:f}  \n').format(gamma)
    lines[5] = ('{:f}  \n').format(gamma)

    # Write the lines back out
    with open('./neboulainput/bc_wall.in', 'w') as f:
        f.writelines(lines)

    #---------------------------------
    #Write input for BL code on out conditions

    with open('./cerbereinput/input.in','r') as f:
        lines = f.readlines()

    # Replace the lines
    lines[1] = ('{:f} \n').format(H) #Hs
    lines[2] = ('{:f} \n').format(Ps) #Ps
    lines[3] = ('{:f} \n').format(Pd) #Pd
    lines[4] = ('{:f}  {:f} \n').format(0.025,Tw) #Reff,Tw

    # Write the lines back out
    with open('./cerbereinput/input.in', 'w') as f:
	    f.writelines(lines)

    # ---------------------------------
    # Execute code to write Neboula input

    call(["../../cerboula/exe/cerbere.exe"])

    # ---------------------------------
    # Execute BL code

    call(["../../cerboula/exe/neboula.exe"])

    # ---------------------------------
    # Read output from BL code

    with open('./neboulaoutput/heatskin.dat') as f:
        lines=f.readlines()
        for line in lines:
            qw_array = np.fromstring(line, dtype=float, sep=' ')

    qw=qw_array[7]

    print(qw)

    return qw

def set_NDPs(NDP1,NDP2,NDP3,NDP4,NDP5):

    #---------------------------------
    # Write gamma for BL code on wall conditions

    with open('./cerbereinput/input.in', 'r') as f:
        lines = f.readlines()

    # Replace the lines
    lines[5] = ('{:f}  {:f}  {:f}  {:f}  {:f}  \n').format(NDP1,NDP2,NDP3,NDP4,NDP5)

    # Write the lines back out
    with open('./cerbereinput/input.in', 'w') as f:
        f.writelines(lines)

    return


def log_likelihood(X): # X = [0,1]

    for i in range(len(X)):
        if X[i]<0. or X[i]>1:
            return -1.e16

    measurements = {"HF": {"mean": 1660000.02, "std-dev": 83000.}, # 10% uncertainty
                    "Pd": {"mean": 33.3, "std-dev": 1.2},
                    "Ps": {"mean": 10000, "std-dev": 0.5},
                    "Tw": {"mean": 350, "std-dev": 17.5},

    }

    priors = {"H": [10.,40.],
              "Pd": [measurements["Pd"]["mean"] - (4*measurements["Pd"]["std-dev"]),measurements["Pd"]["mean"] + (4*measurements["Pd"]["std-dev"])],
              "Ps": [measurements["Ps"]["mean"] - (4*measurements["Ps"]["std-dev"]),measurements["Ps"]["mean"] + (4*measurements["Ps"]["std-dev"])],
              "Tw": [300.,500.],
              "gamma": [-2.,0.]

    }

    H = priors["H"][0] + (X[0]*(priors["H"][1] - priors["H"][0]))
    Tw = priors["Tw"][0] + (X[1]*(priors["Tw"][1] - priors["Tw"][0]))
    Pd = priors["Pd"][0] + (X[2]*(priors["Pd"][1] - priors["Pd"][0]))
    Ps = priors["Ps"][0] + (X[3]*(priors["Ps"][1] - priors["Ps"][0]))
    gamma = priors["gamma"][0] + (X[4]*(priors["gamma"][1] - priors["gamma"][0]))

    gamma = 10**gamma
    H = H*1.e+06

    L = np.divide(np.absolute(measurements["HF"]["mean"] - BL(H,Tw,Pd,Ps,gamma,start_from_previous=False)) 	 ** 2, 2 * measurements["HF"]["std-dev"] ** 2)+np.divide(np.absolute(measurements["Tw"]["mean"]-Tw) ** 2, 2 * measurements["Tw"]["std-dev"] ** 2)+np.divide(np.absolute(measurements["Pd"]["mean"]-Pd) ** 2, 2 * measurements["Pd"]["std-dev"] ** 2)+np.divide(np.absolute(measurements["Ps"]["mean"]-Ps) ** 2, 2 * measurements["Ps"]["std-dev"] ** 2)


    return -L

def m_log_likelihood(X):
    return -1*log_likelihood(X)

NDPs = {"NDP1": 0.215,
        "NDP2": 0.344,
        "NDP3": 0.993,
        "NDP4": 0.213,
        "NDP5": 0.335,

}
set_NDPs(NDPs["NDP1"],NDPs["NDP2"],NDPs["NDP3"],NDPs["NDP4"],NDPs["NDP5"])

## Looking for the MAP point to start sampling ##
n_params = 5
# Xi = [0.5]*n_params
Xi = [0.6, 0.25, 0.32, 0.5, 0.25] # Initial sampling point, choose it based on the Scurve
# res = scipy.optimize.minimize(m_log_likelihood,Xi,method='Nelder-Mead',tol=1e-6)
# print("MAP found at: "+str(res.x))

# MCMC sampling
nburn = 1000
sampler = mcmc.metropolis(np.identity(n_params)*0.01,log_likelihood,nburn)

sampler.seed(Xi)
sampler.Burn()

nchain = 10000
XMCMC = np.zeros((nchain,n_params))

for i in range(nchain):
    XMCMC[i] = sampler.DoStep(1)
    print("MCMC step: "+str(i))

# Saving the chain to external file
filename = './chain_10000.dat'

with open(filename,'w') as ch:
    for i in range(len(XMCMC)):
        for j in range(n_params-1):
            ch.write(str(XMCMC[i,j])+' ')
        ch.write(str(XMCMC[i,n_params-1])+'\n')
