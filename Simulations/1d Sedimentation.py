
# coding: utf-8

# # 1D sedimentation curve
# Using upwind method

# In[4]:

import numpy as np
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
plt.rc('lines',linewidth=0.1)


# #### Grid definition and stability
# I want to set the CFL number to equal 1 to keep stability. I'll change the definition of dt later for more realistic $\Delta t$, $\Delta x$, and $u$.

# In[91]:

# set up time and space steps
nxs = 50 # number of spatial points
L = 10 # height of settling tank
dx = L/(nxs)
xs = np.arange(0+0.5*dx,L+0.5*dx,dx)

c = 5 # m/s wave speed 
cfl = 1 # stability criterion := c*dx/dt
tend = 20 #seconds
dt = dx*c*cfl
nts = int(tend/dt)
ts = np.arange(0,tend,dt)


# #### Setting up initial conditions
# Use a constant concentration profile along the length of the settler
# $$\rho_0 = f(x,0)$$ 
# I need to initialize a storage array too

# In[92]:

u_max = 1
rho_max = 0.95
rho_0 = np.ones_like(xs)*0.5


# ### Define velocity functions

# In[93]:

def fluxfunc(rho,rho_max,u_max,n=1):
    velocity = u_max*(1)
    return rho*velocity

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


# In[95]:

np.shape(ftbs(rho_0,nts,dt,dx,rho_max,u_max))


# In[96]:




# In[ ]:



