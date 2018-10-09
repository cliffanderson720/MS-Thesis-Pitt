# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 13:53:41 2018

@author: cdhig
"""
from velocity_functions import *
import matplotlib.pyplot as plt
import numpy as np

def EOflux(rho):
    fplus = rho*max(RZfluxprime2(rho,1,2.71),0)
    fminus = rho*min(RZfluxprime2(rho,1,2.71),0)
    return fplus,fminus

n=51
rhos = np.linspace(0,1,n)
EOmax  = np.empty(n)
EOmin  = np.empty(n)


for i in range(n):
    EOmax[i]=EOflux(rhos[i])[0]
    EOmin[i]=EOflux(rhos[i])[1]

plt.plot(rhos,EOmin)
plt.plot(rhos,EOmax)

plt.plot(rhos,rhos*RZ(rhos,1,2.71))
