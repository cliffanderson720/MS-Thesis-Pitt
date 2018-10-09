
# coding: utf-8

# In[1]:

import numpy as np
import matplotlib.pyplot as plt


# ### Setting up Variables
# #### Constants:
# > density <br>
# > cell radius <br>
# > stokes dilution coefficients <br>
# > viscosity <br>
# 
# ####  Iterated
# > $C_{total} = \sum{C_i} = 1 $ where i = RBC, WBC, platelets, bacteria <br>
# > $v_i$

# In[2]:

#               g/cm3          cm          seconds           cP
RBC =   {"den": 1.093, "rad": 3.75e-4, "S" : 12e-7,   "visc": False}
WBC =   {"den": 1.066, "rad": 6.25e-4, "S" : 1.2e-7,  "visc": False}
Plate = {"den": 1.053, "rad": 1.19e-4, "S" : 0.032e-7,"visc": False}
Plas  = {"den": 1.024, "rad": False,   "S" : False,   "visc": 2.6}


# #### Operating parameters

# In[3]:

RPM = 3000
omega = RPM/60*2*np.pi # revolutions per second
r = 10 # cm


# #### Variable storage and initialization

# In[4]:

C = 1.0
Psus = 1.024 # arbitrary value

# Time discretization
dt = 0.1                  # time step  (seconds)
tf = 10                   # final time (minutes)
ts = 10                   # number of time points
t  = np.linspace(0,tf*60,ts)

# Spatial discretization
xs = 11                   # number of discrete time points
x  = np.linspace(0,10,xs) # spatial discretization
dx = x[1]-x[0]

# Volume fractions
Crbc = np.empty((ts,xs))  # rows: temporal positions
                          # columns: spacial positions
Cwbc = np.empty((ts,xs))                        
Cplt = np.empty((ts,xs))                        
Cbac = np.empty((ts,xs))
          
# set initial conditions
Crbc[0,:] = 0.4
Cwbc[0,:] = 0.01
Cplt[0,:] = 0.02

Cs = np.dstack((Crbc,Cwbc,Cplt))
for i in range(np.shape(Cs)[2]):
    plt.scatter(x,Cs[0,:,i])


# In[5]:

# test summing along axes
test = np.array([[1,2,3],
                 [2,3,4]])

np.shape(test)[0]


# In[6]:

def vel_i(Ctot,den,rad,S,visc):
    '''
    Calculates the local velocity of the ith cell type
    Inputs: Ctot - previous C_total value(s)
    '''
    
    force   = S*omega**2*r     # centrifugal acceleration term    
    del_rho = (den-Psus)/(den-Plas['den']) # density gradient term
    hinder  = np.exp(-2.5*Ctot/(1-39/64*Ctot)) # local concentration correction introduced by Hawley
    
    return force*del_rho*hinder

vel_rbc = vel_i(Crbc[0,:],**RBC)
vel_rbc


