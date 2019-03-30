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
p1 = np.linspace(0,1)
p2 = 0.00
m = 2
# compute Jacobian?
df1dp1 = p1*(1-p1)**m*(m+1)*(-1) + (1-p1)**(m+1)
df2dp1 = p2*(1-p1)**(m-1)*(-(1/2.)*m+m*p1+p1-1)    #(1+a[1]*m-2*p1)
df2dp2 = (1-p1)**m*(a[1]-p1)


plt.figure(figsize=(10,6))
[plt.plot(p1,dfdx) for dfdx in (df1dp1,df2dp1,df2dp2)]
# plt.plot(p1,df2dp1+df2dp2)
plt.axhline(y=0,color='black',linestyle='--',linewidth=0.5)
# plt.xlim([0,1])
plt.xlabel('Species 1')
plt.ylabel('Jacobian')
plt.legend(['dF1dp1 - row 1 col 1',
            'dF2dp1 - row 2 col 1',
            'dF2dp2 - row 2 col 2'])


# Question: What does the characteristic space look like for single and multi component?
#%% Test full and reduced model for eigenvalues, RI, and waves at all p1, some p2
