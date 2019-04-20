# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 11:06:37 2019

@author: cdhig
"""

#import velocity_functions as vf
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
#%%

def hsc_tot(r1,r2,n=2):
    return (1-r1-r2)**n

def calc_vabs(urel,r1,r2,func=hsc_tot,**kwargs):
    '''return v abs'''
    hsc = hsc_tot(r1,r2,**kwargs)
    urels = np.hstack([1,urel])
    vs  = urels*hsc
    vabs = vs-(np.array([r1,r2])*vs).sum()
    return vabs[-1]

def finda0(r1,r2,func=hsc_tot,**kwargs):
    '''return v abs'''
    hsc = hsc_tot(r1,r2,**kwargs)
    a0 = r1*hsc/(hsc*(1-r2))
    
    return a0

urel = 1/30.

#%%
# define 1 array of r and a plasma volume fraction
nr = 6


aguess = 1/30.
a0 = np.zeros((nr,nr))



#for ax,urel in enumerate(1/np.array([1,2,5,10,15])):
for rp in (0.,0.1,0.25,0.5,0.75):
    r1 = np.linspace(0,1-rp,nr)
    r2 = 1-rp-r1
    a0 = np.zeros()
    
    
    #%%
    vabs = np.empty((nr,nr))
    qabs = np.empty((nr,nr))
    for i in range(nr):
        for j in range(nr):
            vabs[i,j] = calc_vabs(urel,r1[i],r2[j])
            qabs[i,j] = vabs[i,j]*r2[j]

#%%
#    vabs = vabs>=0        
    R1,R2 = np.meshgrid(r2,r1)
#    levels = -np.array([0.1,0.05,0,-.01,-.02,-.05,-.1,-.25,-.5,-.75,-1,-10])
    plt.figure(figsize=(12,4))
    plt.subplot(1,2,1)
    plt.title(f'V_abs for species 2 at rp: {rp:.2f}')
    plt.contourf(R1,R2,vabs,20)
    plt.colorbar()
    
    plt.subplot(1,2,2)
    plt.title(f'Q_abs for species 2 at rp: {rp:.2f}')
    plt.contourf(R1,R2,qabs,20)
    plt.colorbar()

    plt.xlabel('species 1')
    plt.ylabel('species 2')
    plt.show()