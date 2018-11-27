'''
Module containing integration schemes for scalar and system hyperbolic equations.
'''
import numpy as np

def ftbs(rho,nt,dt,dx,rho_max,u_max,flim=None):
    '''
    Forward-time backward-space integration function. 
    It uses the spatial points to the left of a forward-moving wave.
    
    '''   
    rho_n = np.zeros((nt,len(rho)))
    rho_n[0,:] = rho.copy()
    
    for t in range(1,nt):
        F = fluxfunc(u_max,rho_max,rho)
        rho_n[t,1:] = rho[1:] - dt/dx*(F[1:]-F[:-1])
        rho_n[t,0] = rho[0]
        rho = rho_n[t].copy()
    return rho_n  
