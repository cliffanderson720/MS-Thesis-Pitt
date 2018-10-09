"""
Initial crack at integrating the kinematic wave equation
reference: Piao 2017 - platelet recovery paper
"""
#import numpy as np
from numpy import *
#import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import *
from math import *

#u = 1 #m/s
#xs = np.arange(0,5)
#ts = np.arange(0,4)

Ld = 1.0
nTauRun = 0.5
nxy = 22
alpha = 1
cfl = 0.5

tau = Ld**2/alpha
tend = nTauRun*tau
dxy = Ld/(nxy-1)
dt = dxy**2*cfl/(alpha*4)
nt = ceil(tend/dt)
dt = tend/nt
nps = ceil(1/cfl)*10
          
f = zeros((nt,nxy,nxy))
S = ones((nt,nxy,nxy))
rands= random.randint(2,200,size=(50,nxy,nxy))
S[750:800,:,:] = rands

X,Y = meshgrid(linspace(0,Ld,nxy),linspace(0,Ld,nxy))
i = arange(1,nxy-1)
j = i

for it in range(1,nt):
    f[it][ix_(i,j)] = f[it-1][ix_(i,j)]\
    + (alpha*dt/dxy**2)*(f[it-1][ix_(i-1,j)]-2*f[it-1][ix_(i,j)] + f[it-1][ix_(i+1,j)]) \
    + (alpha*dt/dxy**2)*(f[it-1][ix_(i,j-1)]-2*f[it-1][ix_(i,j)] + f[it-1][ix_(i,j+1)]) \
    + S[it-1][ix_(i,j)]
       
    if(it%nps==0):
        clf()
        contourf(X,Y,f[it,:,:],30)
        ion()
        colorbar()
        pause(.001)
        print(it)
#        aa = input()

#contourf(X,Y,f[it,:,:],20)











'''
ran = np.random.random((5,4))
index = np.ix_(xs[::2],ts[::2])
ran
ran[index]
'''