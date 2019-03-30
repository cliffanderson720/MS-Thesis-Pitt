#!/usr/bin/env python
# encoding: utf-8

r"""
1D sedimentation simulation
==================

Solve a one-dimensional polydisperse sedimentation problem

.. math::
    \rho_t + (f(\rho))_x & = 0 \\

Here \rho is a vector of volume fractions \in[0,1]. f is a user-specified function
for the flux of particles.
"""
from __future__ import absolute_import
import importlib
# importlib.reload(riemann)

import numpy as np
from clawpack import riemann
import clawpack as cp

def lowerwall(state,dim,t,qbc,auxbc,num_ghost):
    for i in range(num_ghost):
        qbc[:,i,...] = -qbc[:,2*num_ghost-1-i,...]

def upperwall(state,dim,t,qbc,auxbc,num_ghost):
    for i in range(num_ghost):
        qbc[:,-i-1,...] = -qbc[:,-2*num_ghost+i,...]

def lowerdirichlet(state,dim,t,qbc,auxbc,num_ghost):
    nspecies = qbc.shape[0]
    bvals = np.array([0.,0.])
    for i in range(num_ghost):
        qbc[:,i] = bvals

def upperdirichlet(state,dim,t,qbc,auxbc,num_ghost):
    nspecies = qbc.shape[0]
    bvals = np.array([0.95,0.0])
    for i in range(num_ghost):
        qbc[:,-i-1] = bvals


def setup(use_petsc=False,kernel_language='Python',outdir='./_output',solver_type='classic'):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if kernel_language == 'Python':
        rs = riemann.advection_nonlinear_1D_py.advection_nonlinear_1D
    elif kernel_language == 'Fortran':
        print('No fortran solver available for advection_nonlinear_1D')
        pass

    if solver_type == 'classic':
        solver = pyclaw.ClawSolver1D(rs)
        solver.limiters = pyclaw.limiters.tvd.vanleer
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver1D(rs)

    solver.kernel_language = kernel_language

    solver.bc_lower[0] = pyclaw.BC.custom
    solver.bc_upper[0] = pyclaw.BC.custom
    solver.user_bc_lower = lowerdirichlet
    solver.user_bc_upper = upperdirichlet

    xlower = 0.0
    xupper = 1.0
    mx = 51
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    domain = pyclaw.Domain(x)
    num_eqn = 2
    state = pyclaw.State(domain,num_eqn)

    # Gravitational constant
    state.problem_data['u_rel'] = np.array([1.,1/30.])
    state.problem_data['efix'] = False

    xc = state.grid.x.centers

    IC = 'dam-break'
    # IC = 'uniform-all'
    # IC = 'perturbation'
    x0 = xc[2]

    if IC=='uniform-all':
        c0 = np.array([0.2,0.0])
        # state defaults to empty. Convert to ones and fill with c0
        state.q = np.ones_like(state.q)*c0[:,np.newaxis]

    elif IC=='dam-break':
        # I changed state.is_valid() to always return true for fortran contiguity
        cr0 = np.array([0.2,0.0])
        cl0 = np.array([0.0,0.0])
        state.q = np.ones_like(state.q)
        state.q = cl0[:,np.newaxis]*(xc <= x0)[np.newaxis,:] + \
                  cr0[:,np.newaxis]*(xc >  x0)[np.newaxis,:]
        state.q[0,-1] = 1.

# Change these later to reflect initial conditions
    # elif IC=='2-shock':
    #     hl = 1.
    #     ul = 1.
    #     hr = 1.
    #     ur = -1.
    #     state.q[depth,:] = hl * (xc <= x0) + hr * (xc > x0)
    #     state.q[momentum,:] = hl*ul * (xc <= x0) + hr*ur * (xc > x0)
    elif IC=='perturbation':
        # x1 = x0
        x1 = 0.3
        x2 = 0.7
        eps = 0.1
        state.q[0,:] = eps*np.exp(-1/eps*(xc-x1)**2)
        state.q[1,:] = eps*np.exp(-1/eps*(xc-x1)**2)

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.num_output_times = 50
    claw.tfinal = 10
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.setplot = setplot

    return claw


#-------------------------
def setplot(plotdata):
#--------------------------
    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """
    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for depth
    plotfigure = plotdata.new_plotfigure(name='Component 1', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,1.0]
    plotaxes.ylimits = [0,1.0]
    plotaxes.title = 'Component 1 phi'
    plotaxes.axescmd = 'subplot(211)'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}

    # Figure for species 2[1]
    # plotfigure = plotdata.new_plotfigure(name='Component 2', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = [0,1.0]
    plotaxes.title = 'Component 2'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 1
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}

    return plotdata


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
