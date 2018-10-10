
import numpy as np
import matplotlib.pyplot as plt

def michaels(a,n,amax=0.95):
    return (1-a)**2*(1-a/amax)**n

ns = np.linspace(2,3,6)
a = np.linspace(0,.95)

fig, axs = plt.subplots(2,sharey=True)

for i in range(len(axs)):
    m = 1
    if i==1:
        m = 1-a 
    for n in ns:
        axs[i].plot(a,a*michaels(a,n)/m,label=n)    
    axs[i].legend()
