#!/usr/bin/env python
# encoding: utf-8
r"""
Riemann solvers for a sedimentation problem
"""

from __future__ import absolute_import
import numpy as np
from six.moves import range
import matplotlib.pyplot as plt

# num_eqn = 2
# num_waves = 2

def abs_flux(r,u_rel,hsc=None):
    '''Returns absolute fluxes calculated by Davis 1994 using custom hsc
    Summarizes work in multicomponent_sedimentation.ipynb
    ----------------------------------------------------------------------------
    Inputs: a: array of (j,) stokes velocities relative to first component
            r: array of (i,j) concentrations. Dim (j,) for use in riemann solver
               and in function getJac
            hsc: (optional) default: none. Function for computing hsc. A filler
            function (1-r1-r2)**2 is used.
    ----------------------------------------------------------------------------
    Outputs: qab: (i,j) array of species concentration flux in stationary frame'''
    if r.ndim == 1:
        axis=0
    elif r.ndim > 1:
        axis=1
    else:
        print(f'dimension error in abs_flux. Array r has shape {r.shape}')

    # calculate hindered settling correlation
    hsc = (1-np.sum(r,axis=axis))**2

    if r.ndim > 1: #
        # Do this for (1-r1-r2)**2. make hsc for each species.
        hsc = np.column_stack([hsc for col in range(r.shape[1])])
        cslip = u_rel*hsc
        cab = cslip - np.sum(r*cslip,axis=axis)[:,np.newaxis]
    elif r.ndim == 1:
        cslip = u_rel*hsc
        cab = cslip - np.sum(r*cslip)
        # cab = hsc*(np.ones_like(r)*a-np.sum(a*r,axis=1)[:,np.newaxis]) This one works but sucks
    else:
        print(f'dimenstion error in abs_flux. Array r has shape {r.shape}')
    qab = r*cab
    return qab

def getJac(func,x,*args,**kwargs) :
    """
    Return the numerical Jacobian matrix of a vector x: shape (n,).
    I should maybe allow for x to have ndim > 1 so the function can be called
    without a loop for ndarrays.
    ----------------------------------------------------------------------------
    Inputs:
        func: function that returns dim (n,) array of states (species concentrations)
        x: array dim (n,) of points at which to evalute the jacobian
    ----------------------------------------------------------------------------
    Outputs:
        J: (n,n) Jacobian matrix evaluated at x
    """
    n  = len(x)
    dx = x*1e-8 # optimal dx to balance truncation and roundoff error
    J0 = func(x,*args,**kwargs)[:,np.newaxis] # shape (n,1) to subtract J0 from each column of J1

    J1 = np.zeros((n,n))
    for i in range(n): # loop through columns
        dx_col = np.zeros_like(dx)
        dx_col[i] = dx[i]
        J1[:,i] = func(x+dx_col,*args,**kwargs)
    J  = (J1-J0)/dx

    # for the sake of my BC, I will return zeros and just throw a warning
    nan_vals = np.isnan(J)
    if nan_vals.any():
        # print(f'Nan encountered in getJac for x={x}. Values replaced by 0')
        J[nan_vals] = 0

    return J

# a = np.array([1,1/30.])
# rL = np.array([0.4,0.1])
# rR = np.array([-0.95,-0.01])
# r = np.column_stack([rL,rR])
# r
# rint = np.mean(r,axis=1)
# rint
# J = getJac(abs_flux,rint,a)
# J
#
# s,R = np.linalg.eig(J)
#
# Ri = np.linalg.inv(R)
# s
# Ri

# Js = [getJac(abs_flux,rs,a) for rs in r]
# s,R_raw = np.linalg.eig(Js)
# Ri_raw = np.linalg.inv(R_raw)
# # get the eigenvectors into the right shape for the solver
# R = np.transpose(R_raw,axes=(1,2,0))
# Ri = np.transpose(Ri_raw,axes=(1,2,0))
#
# #%% ---------------------------------------------------------------------------
# problem_data = dict(u_rel = a, efix = False)
# q = np.vstack((r1,r2))

def advection_nonlinear_1D(q_l,q_r,aux_l,aux_r,problem_data):
    r"""
    Polydisperse sedimentation solver in 1d
    no aux arrays are needed
    problem_data may need to contain 'function' key for which hsc to use.
    problem_data will also contain an array of relative velocities.
    """

    # Problem dimensions
    # q_l and q_r must have shapes (num_eqn (or num_waves), num_rp)
    num_rp = q_l.shape[1]
    num_waves = q_l.shape[0]
    num_eqn = q_l.shape[0]
    # I changed the wall BC to reflect all rows. Replaced self.reflect_index[idim] with np.arange(self.num_eqn) in pyclaw/solver._bc_upper

    # Return values
    wave = np.empty( (num_eqn, num_waves, num_rp) )
    s = np.empty( (num_waves, num_rp) )
    amdq = np.zeros( (num_eqn, num_rp) )
    apdq = np.zeros( (num_eqn, num_rp) )

    # Solver parameters
    u_rel = problem_data['u_rel']

    # Calculate average states
    # This is not the Roe average. Problems may occur around shocks or boundaries
    q_ave = (0.5*q_l+0.5*q_r)
    # q_ave = q_r
    # manually set dirichlet condition on first upper ghost cells after avg.
    q_ave[:,-2] = q_ave[:,-1]
    # manually set wall on first upper ghost cell
    # q_ave[:,-2] = -q_ave[:,-3]
    # q_ave[:,-1] = -q_ave[:,-4]

    if (q_ave>1).any():
        print('value in q_ave > 1')
        q_ave[q_ave>1]=1.


    # Evaluate jacobian - when checking this for consistency, remember q_ave != q_l
    Jac = [getJac(abs_flux,q_ave[:,i],u_rel) for i in range(num_rp)]

    # Find eigenvector coefficients of jacobian at each interface
    s_raw,R_raw = np.linalg.eig(Jac)
    Ri_raw      = np.linalg.inv(R_raw)

    # reshape the R arrays to match the shape of the wave variable
    # slicing R as [:,:,i] gives the eigenvector matrix of the ith interface
    R  = np.transpose(R_raw,axes=(1,2,0))
    Ri = np.transpose(Ri_raw,axes=(1,2,0))
    s  = s_raw.T

    # test eigevvectors for complex Values
    if np.iscomplexobj(R):
        print('Complex value encountered in eigenvectors')

    # Get riemann invariants (Rinv*(q_r-q_l))
    delta = q_r - q_l

    # removed hardcoding the wave variable
    # Compute the waves
    wave = Ri*delta
    # np.set_printoptions(precision=3,linewidth=np.inf)
    # Entropy fix
    if problem_data['efix']:
        raise NotImplementedError("Entropy fix does not exist for solver advection_nonlinear_1D!")
    else:
        # Godunov update
        s_index = np.zeros((2,num_rp))
        for m in range(num_eqn):
            for mw in range(num_waves):
                s_index[0,:] = s[mw,:]
                # There may still be a broadcasting error or something,
                # but if there is, I can't trace it.
                # print(m,mw)
                # print(f'wave[{m,mw}]')
                # print(wave[m,mw,:])
                # print(f's[{mw},:]')
                # print(s[mw,:])
                # print('min (s_index,0)')
                # print(np.min(s_index,axis=0))
                # print('product smin*wave')
                # print(np.min(s_index,axis=0) * wave[m,mw,:])

                amdq[m,:] += np.min(s_index,axis=0) * wave[m,mw,:]
                apdq[m,:] += np.max(s_index,axis=0) * wave[m,mw,:]
                # print(amdq)
                # print(apdq)
                # print()
    # amdq[:,11:]=0.
    # np.set_printoptions(linewidth=np.inf)
    # print('wave speed')
    # print(s)
    #
    # print('ql,qr,qave')
    # print(q_l)
    # print(q_r)
    print('q_ave at interface')
    print(q_ave)
    # print('waves')
    # print(q_ave[:,8:])
    # [print(i,'\n',wave[...,i],s[:,i],'\n',np.column_stack([apdq[:,i],amdq[:,i]])) for i in (8,9,10,11,12,13)]
    print('amdq')
    print(amdq)
    print('apdq')
    print(apdq)
    return wave,s,amdq,apdq


# advection_nonlinear_1D(q[:,1:],q[:,:-1],None,None,problem_data)
# fig,ax = plt.subplots(nrows=2,ncols=2)
