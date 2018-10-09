# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 10:37:22 2018

@author: cdhig
"""
import numpy as np
import matplotlib.pyplot as plt


def volflux(a,k=0.65,m=2.5):
    return a*(1-a)**2/((1-a/k)**(-k*m))

a = np.linspace(0,0.7)

plt.plot(a,volflux(a),label='k=0.65,m=2.5')
plt.plot(a,volflux(a,k=1,m=3),label='k=1,m=3')
plt.legend()
