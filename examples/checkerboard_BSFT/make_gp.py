#This code is part of DIETA
#
#Authored by Kirill Alpin

import numpy as np
from dieta_wrapper import DIETA

#define the chemical potential and hopping amplitude
mu = 1.0
t = -0.5

#define the parameters given to DIETA
params = {
    "mu0" : mu,
    "mu1" : -mu,
    "mup0" : 0.0,
    "mup1" : 0.0,
    "t" : t
    }

#create an instance of dieta with 
#the graph file found in graph_2site.txt
dieta = DIETA("graph_2site.txt", params)


#plotting part
import matplotlib.pyplot as plt

D=[]
Omega=[]
#iterate over Delta
for Delta in np.arange(-mu, mu/2.0, 0.01):
    D.append(Delta)
    
    #vary the chemical potentials
    params["mu0"] = mu + Delta
    params["mu1"] = -mu + Delta
    params["mup0"] = -Delta
    params["mup1"] = -Delta
    
    #update the parameters
    dieta.update_params(params)
    
    #calculate the BSFT functional
    Omega.append(dieta.get_BSFT_functional())
    
#plot the curve D, Omega
lt='-'
plt.plot(D, Omega, lt)
plt.xlabel("$\Delta$")
plt.ylabel("$\Omega(\Delta)$")
plt.show()
plt.savefig('gp_plot.pdf')
