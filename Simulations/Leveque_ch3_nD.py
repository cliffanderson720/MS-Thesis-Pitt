# -*- coding: utf-8 -*-
"""
Created on Sat Feb  2 13:32:01 2019

@author: cdhig
"""

import numpy as np
import matplotlib.pyplot as plt

def abs_flux(a,r,hsc=None):
    '''Returns absolute fluxes calculated by Davis 1994 using custom hsc
    Summarizes work in multicomponent_sedimentation.ipynb
    Inputs: a: array of (j,) stokes velocities relative to first component
            r: array of (i,j) concentrations
            hsc: (optional) default: none. Function for computing hsc. A default
            function (1-r1-r2)**2 is used.
    Outputs: qab: (i,j) array of species concentration flux in stationary frame
    THIS CODE LIVES IN THE advection_nonlinear_1D_py.py FILE IN THE RIEMANN PACKAGE'''
    # if r.ndim == 1: # if 1D array, sum
    #     axis=0
    # elif r.ndim>1:
    #     axis=1
    # else:
    #     print(f'dimension problem in function {__name__}')

    axis=1
    hsc1d = (1-np.sum(r,axis=axis))**2
    hsc = np.column_stack([hsc1d for col in range(r.shape[1])])
    # broadcasting nightmare. Simplify this, please
    cslip = a*hsc
    cab = cslip - np.sum(r*cslip,axis=axis)[:,np.newaxis]
    # cab = hsc*(np.ones_like(r)*a-np.sum(a*r,axis=1)[:,np.newaxis]) This one works but sucks
    qab = r*cab
    return qab


r1 = np.linspace(0,0.95,6)
r2 = np.linspace(0,0.05,6)
r = np.column_stack((r1,r2))
a = np.array([1,1/30.])
fluxes = abs_flux(a,r)
fluxes

def getJac(func,x) :
    """
    Return the numerical Jacobian matrix
    Inputs:
        func: function that outputs a nxm vector of function values
        x:    1D array of values where the Jacobian will be evaluated
    Outputs:
        J: values for the Jacobian (n,n)

    """
    n  = len(x)
    dx = x*1e-8
    J0 = func(x)[:,np.newaxis]  # shape (n,1) to subtract J0 from each column of J1

    J1 = np.zeros((n,n))
    for i in range(n): # loop through columns
        dx_col = np.zeros_like(dx)
        dx_col[i] = dx[i]
        J1[:,i] = func(x+dx_col)

    J  = (J1-J0)/dx

    return J

getJac(fluxes)





# A = np.array([[-4,3],
#               [1,2]])
# ql = np.array([1,3])
# qr = np.array([0,-1])



#A = np.array([[-4,3,1],
#              [1,2,1],
#              [7,-1,0]])
#ql = np.array([1,3,1])
#qr = np.array([0,-1,-7])
# Eigenvalue decomposition
lamda,R = np.linalg.eig(A)
L  = np.diag(lamda)
Ri = np.linalg.inv(R)

# get characteristic values on left and right side
wl = np.dot(Ri,ql)
wr = np.dot(Ri,qr)

alpha = wr-wl
W = alpha[np.newaxis,:]*R # multiply columns of R by alpha
# Get columns of W where lamda <0
# use np.where()[0] toretrieve 1D array and broadcast to column numbers
Wneg = W[:,np.where(lamda<0)[0]]
qstar = ql+Wneg.sum(axis=1)

Wpos = W[:,np.where(lamda>=0)[0]]
qstar = qr-Wpos.sum(axis=1)


def constructq(w,R):
    ''' w is an array of wp values shape: (n,)
        R is a matrix of right eigenvectors shape: (n,n)'''
    # broadcast w to have shape: (1,n)
    # multiply columns of R by new w and sum the columns
    return np.sum(w[np.newaxis,:]*R,axis=1)
#%% Make graphs
# make graphs of characteristics
x  = np.arange(2)
[plt.plot(x*lam,x) for lam in lamda]
plt.xlabel('x')
plt.ylabel('t')
plt.show()

#%%
# Make graphs of hugoniot loci
# plot the value of ql and qr in the (q1,q2) plane. Draw eigenvectors
# and find qstar
qlr = np.vstack((ql,qr)).T
plt.plot(qlr[0,:],qlr[1,:],'o')
# draw a line between 0 and rp
z = np.zeros(2)
r1 = np.vstack((z,R[:,0]))
r2 = np.vstack((z,R[:,1]))

plt.plot(r1[:,0],r1[:,1],'k') # unit r1
plt.plot(r2[:,0],r2[:,1],'k--') # unit r2
# ql = wl1*r1 + wl2*r2
ql_fig = constructq(wl,R)
qr_fig = constructq(wr,R)
qlr_fig = np.vstack((ql_fig,qr_fig)).T
plt.plot(qlr_fig[0,:],qlr_fig[1,:],'bo')

# plot the intermediate state
plt.plot(qstar[0],qstar[1],'go')
plt.ylabel('q2')
plt.xlabel('q1')
