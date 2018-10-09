#src-ch6/advection_schemes.py
#http://folk.ntnu.no/leifh/teaching/tkt4140/._main072.html#ch:6_Lax-Wendroff_nonlin
# I took this code from the URL above because it had a lot of the structure I wanted

import numpy as np
from matplotlib import animation
from scipy import interpolate
import matplotlib.pylab as plt
from velocity_functions import *

radius  = {'RBC':3.75e-6,'WBC':6.25e-6,'plate':1.19e-6,'bac':2.5e-6} # radius in meters
density = {'RBC':1093,   'WBC':1066,   'plate':1053,   'bac': 1025 } # density in kg/m3
sedcoef = {'RBC':12e-7,  'WBC':1.2e-7, 'plate':0.032e-7} # Sedimentation coefficient in s (Van Wie)
viscosity = {'plas':0.0026} # viscosity of plasma in kg/m s

init_func=0   # Select stair case function (0) or sin^2 function (1)

# function defining the initial condition
if (init_func==0):
    def f(x):
        """Assigning a value of 1.0 for values less than 0.1"""
        f = np.zeros_like(x)
        f[np.where(x <= 0.042)] = 1.0
        f[np.where(x >  0.042)] = 0.1
        return f
elif(init_func==1):
    def f(x):
        """A smooth sin^2 function between x_left and x_right"""
        f = np.zeros_like(x)
        x_left = 0.25
        x_right = 0.75
        xm = (x_right-x_left)/2.0
        f = where((x>x_left) & (x<x_right), np.sin(np.pi*(x-x_left)/(x_right-x_left))**4,f) 
        return f

def ftbs(u): # forward time backward space
    u[1:-1] = u[1:-1] + c*(u[:-2]-u[1:-1])
    return u[1:-1]

# Lax-Wendroff
def lax_wendroff(u): 
    u[1:-1] = c/2.0*(1+c)*u[:-2] + (1-c**2)*u[1:-1] - c/2.0*(1-c)*u[2:]
    return u[1:-1]

# Lax-Friedrich Flux formulation
def lax_friedrich_Flux(u):
    u[1:-1] = (u[:-2] +u[2:])/2.0 -  dt*(F(u[2:])-F(u[:-2]))/(2.0*dx)
    return u[1:-1] 

# Lax-Friedrich Advection
def lax_friedrich(u):
    u[1:-1] = (u[:-2] +u[2:])/2.0 -  c*(u[2:] - u[:-2])/2.0
    return u[1:-1] 

# macCormack for advection quation
def macCormack(u):
    up = u.copy()
    up[:-1] = u[:-1] - c*(u[1:]-u[:-1])
    u[1:] = .5*(u[1:]+up[1:] -  c*(up[1:]-up[:-1]))
    return u[1:-1] 


# Constants and parameters
amax = stokes(3000,0.06,radius['RBC'],density['RBC']) # wave speed
tmin, tmax = 0.0, 40  # seconds. start and stop time of simulation
xmin, xmax = 0.04, 0.06 # meters. start and end of spatial domain
Nx = 80 # number of spatial points
c = 0.9 # courant number, need c<=1 for stability


# Discretize
x = np.linspace(xmin, xmax, Nx+1) # discretization of space
dx = float((xmax-xmin)/Nx) # spatial step size
dt = c/amax*dx # stable time step calculated from stability requirement
Nt = int((tmax-tmin)/dt) # number of time steps
time = np.linspace(tmin, tmax, Nt) # discretization of time

# solve from tmin to tmax
solvers = [ftbs,lax_wendroff,lax_friedrich,macCormack]
u_solutions = np.zeros((len(time),len(x),len(solvers)))
uanalytical = np.zeros((len(time), len(x))) # holds the analytical solution
    
for k, solver in enumerate(solvers): # Solve for all solvers in list
    u = f(x) # Initial conditions setup. For multicomponent...
    un = np.zeros((len(time), len(x))) # holds the numerical solution

    for i, t in enumerate(time[1:]):
        
        if k==0:
            uanalytical[i,:] = f(x-amax*t) # compute analytical solution for this time step
            
        u_bc = interpolate.interp1d(x[-2:], u[-2:]) # interplate at right bndry
        
        u[1:-1] = solver(u[:]) # calculate numerical solution of interior
        u[-1] = u_bc(x[-1] - amax*dt) # interpolate along a characteristic to find the boundary value
        
        un[i,:] = u[:] # storing the solution for plotting
    
    u_solutions[:,:,k] = un

#%%
X,Time = np.meshgrid(x,time)
plt.contourf(X,Time,u_solutions[:,:,0])
plt.colorbar()

#%%
### Animation 
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(xmin,xmax), ylim=(np.min(un), np.max(un)*1.1))

lines=[]     # list for plot lines for solvers and analytical solutions
legends=[]   # list for legends for solvers and analytical solutions

for solver in solvers:
    line, = ax.plot([], [])
    lines.append(line)
    legends.append(solver.__name__)

line, = ax.plot([], []) #add extra plot line for analytical solution
lines.append(line)
legends.append('Analytical')
time_text = ax.text(0.05,1.05,'')

plt.xlabel('Radius (m)')
plt.ylabel('Volume fraction')
plt.legend(legends, loc=3, frameon=False)


# initialization function: plot the background of each frame
def init():
    for line in lines:
        line.set_data([], [])
    return lines,

# animation function.  This is called sequentially
def animate_alt(i):
    for k, line in enumerate(lines):
        if (k==len(lines)-1):
            line.set_data(x, uanalytical[i,:])
        else:
            line.set_data(x, u_solutions[i,:,k])
            time_text.set_text('t={0:.1f} s'.format(time[i]))
    return lines,time_text

 
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate_alt, init_func=init, frames=Nt, interval=150, blit=False)
plt.show()