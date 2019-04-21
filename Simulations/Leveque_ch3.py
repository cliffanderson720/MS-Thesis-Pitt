# -*- coding: utf-8 -*-
"""
Created on Sat Feb  2 13:32:01 2019

@author: cdhig
"""

import numpy as np
import matplotlib.pyplot as plt

A = np.array([[-4,3],
              [1,2]])
ql = np.array([1,3])
qr = np.array([0,-1])
# Eigenvalue decomposition
lamda,R = np.linalg.eig(A)
L  = np.diag(lamda)
Ri = np.linalg.inv(R)

# get characteristic values on left and right side
wl = np.dot(Ri,ql)
wr = np.dot(Ri,qr)

alpha = wr-wl
W = alpha[np.newaxis,:]*R
# Get columns of W where lamda <0
Wneg = W.T[np.where(lamda<0)]
qstar = ql+Wneg[0]

Wpos = W.T[np.where(lamda>=0)]
qstar = qr-Wpos[0]


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



