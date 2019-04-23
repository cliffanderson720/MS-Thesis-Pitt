# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 14:19:52 2019

@author: cdhig
"""

import numpy as np
A = np.arange(9).reshape((3,3))
A[-1,-1] = 2.
L,R = np.linalg.eig(A)
Rm  = np.linalg.inv(R)
A_recon = np.dot(np.dot(R,np.diag(L)),Rm)

q = np.array([2,3,10])
x = np.linalg.solve(A,q)

#test. Is R Rm q=q the same as q=sum_{p=1}{wp rp}?
W = np.matmul(Rm,q)
W