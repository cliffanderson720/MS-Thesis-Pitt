'''
Module containing the integrator class
'''
import velocity_functions as vf
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

class driver:
    '''
    Sets up and solves 1D scalar sedimentation problem
    '''
    def __init__(self,H0,ngrd=100,L=1,tend=2,ntimes=20,
                 hsc=vf.michaels,cfl=0.95,H_loc = None, a_loc=0.5,a_right=0.99):
        '''Take initial conditions and desired hindered settling correction.
        Sets up computational grid given user inputs and integrates.'''
        # get initial condition and hindered settling function
        self.H0 = H0
        self.hsc = hsc
        
        # set spatial grid-----------------
        self.ngrd  = ngrd
        self.L     = L
        self.ngrdf = self.ngrd+1
        self.Δx    = self.L/self.ngrd
        self.xf    = np.linspace(0,self.L,self.ngrdf)
        self.x     = (self.xf[1:] + self.xf[:-1])/2
        self.a_center = np.ones(self.ngrd)
        self.a_center[int(a_loc*self.ngrd):] = a_right
        self.a_left = 1
        self.a_right = a_right

        ## set time grid--------------------
        self.tend = tend

        # get max dimensionless speed
        self.v0     = self.get_v(self.H0)
        
        # adjust dt if cfl is specified
        # if user has changed CFL, set ntimes using CFL and dx
        if cfl not in self.__init__.__defaults__: 
            self.cfl    = cfl
            self.dt     = self.cfl*self.Δx/self.v0
            self.times  = np.arange(0,self.tend+self.dt,self.dt)
            self.ntimes = len(self.times)

        else:
            self.ntimes = ntimes # else, just use the default ntimes
            self.times = np.linspace(0,self.tend,self.ntimes)
            self.cfl = self.v0*(self.times[1]-self.times[0])/self.Δx
        
        self.dt = self.times[1]-self.times[0]
        
        # Give initial concentration profile and integrate
        self.ϕ0   = np.ones(ngrd)*self.H0
        if H_loc is not None:
            self.ϕ0[int(H_loc*self.ngrd):] = 0
        self.soln = abs(self.integrate())
        return

    def get_v(self,ϕ,volavg=True):
        '''
        Returns velocity for a given ϕ. 
        volavg (Boolean) specifies absolute (True)
        or slip velocity (False).
        '''
        c = np.ones_like(ϕ)
        if volavg:
            v = self.hsc(ϕ)*(1-ϕ)
        else:
            v = self.hsc(ϕ)
        return v

    def get_bc(self,center,left,right):
        '''
        Returns array of values on all grid boundaries  
        '''
        a_bc = np.hstack([left,center,right])
        return a_bc

    def godunov(self,ϕ,t,ϕ_max):
        '''
        Godunov flux function with boundary conditions.
        Sets up dϕdt with entropy-satisfying fluxes in a
        method of lines formulation. 
        Called by method "integrate"
        '''
        # set wall boundary conditions
        ϕ_bc = self.get_bc(ϕ,0,ϕ_max)
        
        # compute flux for each cell using volume-averaged ϕ
        q = self.get_v(ϕ_bc) * ϕ_bc #* self.get_bc(self.a_center,self.a_left,self.a_right)

        # evaluate fluxes at each cell wall 
        qf = np.zeros(self.ngrdf)
        for i in range(self.ngrdf):
            if ϕ_bc[i] >= ϕ_bc[i+1]:
                qf[i] = max(q[i],q[i+1])
            elif ϕ_bc[i] <= 0 and 0 <= ϕ_bc[i+1]:
                qf[i] = 0
            else:
                qf[i] = min(q[i],q[i+1])
        
        # Method of lines form
        dϕdt = 1/self.Δx*(qf[:-1]-qf[1:])
        return dϕdt
    
    def plot(self,plot_type='snap'):
        '''
        Plots snapshots of the concentration profile with time. 
        Can plot snapshots in time or a gif
        '''
        if plot_type=='snap':
            color  = plt.get_cmap('Reds')
            fig,ax = plt.subplots(ncols=1)
            colors = iter(color(np.linspace(0,1,self.ntimes)))
            for j in range(self.ntimes):
                ax.plot(self.x,self.soln[j,:],color=next(colors))
            plt.show()
        
        elif plot_type=='gif':
            from matplotlib import animation
            plt.rc('animation',html='html5')

            fig,ax = plt.subplots()
            
            ax.set_xlim(left=0, right=1)
            ax.set_ylim(bottom=-0,top=1.0)
            ax.set_xlabel(r'Dimensionless position $r^\ast$')
            ax.set_ylabel(r'Hematocrit ($\phi$)')
            
            line, = ax.plot(self.x, self.soln[0,:], lw=2,color='r')
            time  = ax.text(0.1,0.5,'hi')
            
            def animate(i):
                line.set_data(self.x,self.soln[i,:])
                time.set_text(f't* = {self.times[i]:0.1f}')
                return (line,time)
            
            anim = animation.FuncAnimation(fig, animate, frames=self.ntimes, interval=self.dt*1000*4.7)
            anim.save('sim.html') #,writer='pillow')
            plt.draw()
            plt.show()
        else:
            raise Exception('plot_type must be a string with value "gif" or "snap"')
        return 

    def integrate(self):
        '''Returns concentration profile at requested time intervals'''
        self.ϕall = odeint(self.godunov, self.ϕ0, self.times,args=(0.95,))
        return self.ϕall
