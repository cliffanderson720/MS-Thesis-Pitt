"""
Created on Mon Jul 23 15:34:23 2018

@author: cdhig
Contains hindered settling functions for sedimentation
"""
import numpy as np
import matplotlib.pyplot as plt

def without_keys(dictionary, keys):
    '''
    Excludes keys from a dictionary
    Inputs:
        dictionary: dictionary object to select from
        keys: tuple of key strings to include

    '''
    return {key: value for key, value in dictionary.items() if key not in keys}

def with_keys(dictionary, keys):
    """
    Filters a dict by only including certain keys
    Inputs:
        dictionary: dictionary object to select from
        keys: tuple of key strings to include
    """
    key_set = set(keys) & set(dictionary.keys())
    return {key: dictionary[key] for key in key_set}

#def select_cells(cell_types,dictionaries = (radius,density,sedcoef,viscosity),action='include'):
#    '''
#    Subsets dictionaries to include desired cell types. This function exists to
#    preclude individual functions from subsetting. This may speed things up.
#    --------------------------------------------------------------------------
#    Inputs:
#        cell_types: tuple of cell types strings to include/exclude. Default action is to include
#                    possible cell types include 'RBC', 'WBC', 'plate', 'bac'
#        dictionaries: tuple of dictionaries to subset
#        action: string selecting whether to 'include' or 'exclude' cell_types
#    --------------------------------------------------------------------------
#    Output:
#        Tuple of dictionaries to use for calculations
#    '''
#    
#    
#    '''
#    I stopped b/c probably not worth it. But here's a good url:
#        https://stackoverflow.com/questions/17659580/python-list-of-dictionaries-projection-filter-or-subset
#    '''
#    if action == 'include':
#        selected_cells = with_keys(dictionary,cell_types)
#    return

radius  = {'RBC':3.75e-6,'WBC':6.25e-6,'plate':1.19e-6,'bac':2.5e-6} # radius in meters
density = {'RBC':1093,   'WBC':1066,   'plate':1053,   'bac': 1025 } # density in kg/m3
sedcoef = {'RBC':12e-7,  'WBC':1.2e-7, 'plate':0.032e-7} # Sedimentation coefficient in s (Van Wie)
viscosity = {'plas':0.0026} # viscosity of plasma in kg/m s

    
def stokesdict(RPM,r,radius_dict,density_dict,rho_f=1024,visc=.0026,mode=None):
    '''
    Calculates the stokes velocity in m/s for j particle types
    
    -------------------------------------------------------------------------
    Inputs: 
        RPM: scalar: revoluctions per minute
        r: scalar: moment arm from the axis of rotation to the cell location (m)
        r_dict: dict: keyword from a dictionary of cell radii. Can be dictionary or key.
        rho_cell: dict: same as r_dict, but density
    Optional:
        mode: string: default None, 'gravity'
    Output:
        sedimentation velocity in m/s
    '''
    # Check to see if the inputs are dictionaries. If dictionaries, unpack
    # if floats, move on
    if all((type(radius_dict)==dict,type(density_dict)==dict)) == True:
        r_cell = np.array(list(radius_dict.values()))
        rho_i = np.array(list(density_dict.values()))
#        print('this is a dictionary')
    else:
        r_cell = radius_dict
        rho_i = density_dict
#        print('this is a float')


    if mode == 'gravity':
        g = 9.8 # m/s2
        return (2*(rho_i-rho_f)*g*r_cell**2)/(9*visc)
    elif mode == 'hawkseley':
        print('hi, this is hawkseley speaking. I am not here right now')
    else:
        omega = RPM*2*3.1416/60  # convert from revolutions per minute to radians/s
        accel = r*omega**2
        return accel*(2*(rho_i-rho_f)*r_cell**2)/(9*visc)


def stokes(RPM,r,cell_radius,cell_density,rho_f=1024,visc=.0026,mode=None, **kwargs):
    '''
    Calculates the stokes velocity in m/s for j particle types
    -------------------------------------------------------------------------
    Required inputs: 
        RPM: scalar: revolutions per minute
        r: ndarray length i: moment arm from the axis of rotation to the cell location (m)
        r_dict: tuple of j cell radii (m)
        rho_cell: tuple of j cell densities (kg/m3)
    Optional inputs:
        mode: string: default None
        Options: 'gravity' - currently not working
    Output:
        speed: ndarray (jxi): stokes velocity of j particle types at i spatial points
    '''
    omega = RPM*2*3.1415926353/60  # convert from revolutions per minute to radians/s
    accel = r*omega**2

    try:
        ncells = len(cell_radius)
        nspatial = len(r)
        speeds = np.empty((ncells,nspatial))
        for j in range(ncells):
            speeds[j] = accel*(2*(cell_density[j]-rho_f)*cell_radius[j]**2)/(9*visc)

    except TypeError:
        speeds = accel*(2*(cell_density-rho_f)*cell_radius**2)/(9*visc)
        
    return speeds



def RZ(rho,umax,n,rhomax=1):
    return umax*(rhomax-rho)**n

def RZfluxprime(rho,umax,n,rhomax=1):
    return umax*((1-n)*(1-rho)**n + n*(1-rho)**(n-1))

def RZfluxprime2(rho,umax,n):
    '''
    Computes the derivative of the RZ volume flux. This expression was derived
    in mid july 2018 by hand by Clifton Anderson.
    -------------------------------------------------------------------------
    inputs: 
        rho: n by j array of particle concentrations
             n: rows: if 1D, n=1. If 2D, n is the length of the input array
             j: cols: number of cell types        
        umax: scalar: stokes velocity for each particle type
        n: scalar >0. Same value as used for RZ flux.
    -------------------------------------------------------------------------
    outputs:
        n by j array of derivative values
    '''
    return -rho*umax*n*(1-rho)**(n-1) + umax*(1-rho)**n

def porosity(concs,power=1,radius_dict=radius):
    '''
    Calculates the apparent porosity for j particle types
    -------------------------------------------------------------------------
    inputs:
        concs: n by j array of particle concentrations. Values are summed by row
               n: rows: if 1D, n=1. If 2D, n is the length of the input array
               j: columns: number of cell types. 

        power: scalar exponent for the porosity. Used in empirical correlations.
               see Masliyah articles quoted in Patwardhan and Tien 1985
        
        radius_dict: dictionary object containing the radius of included cell types
    -------------------------------------------------------------------------
    outputs:
        porosity: n by j array of apparent porosities for each input composition
        
    Notes: When plotting a single row of concentrations, change the shape of
           conc[i,:] from (j,) to (j,1) using np.atleast_2d() or np.reshape()
    '''
    np.seterr(divide='ignore',invalid='ignore')
    
    # convert dictionary of cell values to an array to allow for array operations
    d_cells = 2*np.array(list(radius_dict.values()))

    # check to see if input concentration array is 1D or 2D. If 2D, sum row-wise
    try:
        if len(concs.shape) == 1:
            dim = 0
        elif len(concs.shape) == 2:
            dim = 1
    except:
        print('func: porosity: dimension error in concentration summing')
        
    
#    # normalize concentrations to one
#    concs = concs / np.sum(concs,axis=dim,keepdims=True) # for keepdims explanation, see https://stackoverflow.com/questions/16202348/numpy-divide-row-by-row-sum
    
    # Compute the void envelope for the particle mixture
    d_average = np.sum(concs*d_cells,axis=dim)/np.sum(concs,axis=dim)
    void_envelope = d_average*(np.sum(concs,axis=dim)**(-1/3.)-1)
    
    # reshape the data into n rows and 1 column
    void_envelope = np.reshape(void_envelope,(len(void_envelope),1))
    
    # compute the apparent porosity for each particle type
    # I put the nan function there to make the porosity = 1 for particle concentration = 0
    porosity = 1-np.nan_to_num((1+void_envelope/d_cells)**-3)
    return porosity**power


#### Setting up concentration arrays to test the porosity function
#conc1a = np.linspace(0,.25,51)
#concs = np.empty((4,len(conc1a)))
#concs = np.array([conc1a for i in range(4)]).T
#                
##porosity(concs)
#plt.plot(concs,porosity(concs,power=2.71))
#plt.legend(radius.keys())
####

#RPMs = np.linspace(0,8000)
#fig,axs = plt.subplots(2)
#axs[0].plot(RPMs,stokes(RPMs,0.06,radius['RBC'],density['RBC']))
