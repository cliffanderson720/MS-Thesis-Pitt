'''Checking numerical jacobian values against analytical ones for the
reduced model where all speed is a function of RBC only. Numerical jacobian
results are not recorded, but the call signature is
J = getJac(abs_flux,rL,a,basis=0)
The signs of each element in jacobian depends on the phi 1, a_n1, and phi n
'''
import numpy as np
import matplotlib.pyplot as plt
# testing analytical expression for jacobian on march 5 2019
a = np.array([1,1/30.])
x1 = np.linspace(0,1)
x2 = 0.1
m = 2
df1dx1 = x1*(1-x1)**m*(m+1)*(-1) + (1-x1)**(m+1)
df1dx2 = x2*(1-x1)**(m-1)*(-(1/2.)*m+m*x1+x1-1)    #(1+a[1]*m-2*x1)
df1dx2[0]
-x2*(1+a[1]*m)
df2dx2 = (1-x1)**m*(a[1]-x1)
(df2dx2<0).any()
[plt.plot(x1,dfdx) for dfdx in (df1dx1,df1dx2,df2dx2)]
plt.xlim([0,1])
plt.legend(['dF1dx1 - row 1 col 1',
            'dF1dx2 - row 1 col 2',
            'dF2dx2 - row 2 col 2'])
plt.plot()
