import numpy as np
import matplotlib.pyplot as plt

def moody(a,n,amax=0.95):
    return (1-a)**2/((1-a/amax)**-n)

ns = np.linspace(0,2,6)
a = np.linspace(0,.95)

for n in ns:
    plt.plot(a,moody(a,n),label=n)
    
plt.legend()
plt.xlabel('concentration')