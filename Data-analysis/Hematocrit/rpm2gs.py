import numpy as np
import matplotlib.pyplot as plt

n=201
rpms = np.linspace(0,3500,n)
rs = np.linspace(0,0.15,n) # radius of the rotor in meters

def gfunc(rpm,r):
    accel = rpm*2*np.pi/60
    return r*accel**2

gfunc(rpms,rs)
