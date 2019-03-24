'''
symbolic manipulation to do eigenvalue decomposition on jacobian
'''
from sympy import *
#%%

init_printing(False)
a1,a2,a3,p1,p2,p3,n,hsc = symbols('a1,a2,a3,p1,p2,p3,n,hsc')
p = Matrix([p1,p2]).as_mutable()
p
hsc = (1-p1)**n

flux_abs = Matrix([hsc*p1*(1 -p1),
                   hsc*p2*(a2-p1)]) # for 3rd component: hsc*p3*(a3-p1)])

jac_abs = flux_abs.jacobian(p)
jac_abs

# convert to eigenspace
eispace_raw = jac_abs.eigenvects()
eispace = eispace_raw.copy()

eispace = [item for component in eispace for item in component]
eivals_abs = eispace[0::3]
multiplicity_abs = eispace[1::3]
eivecs_abs = eispace[2::3]

[collect(eival_abs,(-p1+1)**n) for eival_abs in eivals_abs]
eivecs_abs

R = Matrix([[eivecs_abs[i][0][j] for i in range(len(eivecs_abs))] for j in range(len(eivecs_abs))]).as_mutable()
Ri = eivecMat.inverse_GE().as_mutable()
Ri*p #<--plot this as f(p1,p2)

Rfunc = lambdify((p1,p2,a2,n),R,'numpy')

#%% plot eigenvalues
import numpy as np
import matplotlib.pyplot as plt
a2 = 1/30.
nr = 101
r1 = np.linspace(0,0.9,nr)
r2 = 0.1
m = 2
# eigenvalues
lam1 = ((1-r1)**m)*(a2-r1)
lam2 = ((1-r1)**m)*(-m*r1-2*r1+1)

# Riemann invariants:

RI = np.array([Rfunc(r1i,r2,a2,m) for r1i in r1])
RI[0,...]

# waves computed from Ri*p
W1 = r1*(-m*r2*(a2 - r1)*(-r1 + 1)**m/(-r1 + 1) - r2*(-r1 + 1)**m)/(m*r1*(-r1 + 1)**m + r1*(-r1 + 1)**m + (a2 - r1)*(-r1 + 1)**m - (-r1 + 1)*(-r1 + 1)**m) + r2
W2 =-r1*(-m*r2*(a2 - r1)*(-r1 + 1)**m/(-r1 + 1) - r2*(-r1 + 1)**m)/(m*r1*(-r1 + 1)**m + r1*(-r1 + 1)**m + (a2 - r1)*(-r1 + 1)**m - (-r1 + 1)*(-r1 + 1)**m)
# W1[np.where(W1>1)]=1
# W2[np.where(W2<-1)]=-1
W1+W2

fig,ax = plt.subplots(nrows=3,sharex=True)
ax[0].plot(r1,lam1)
ax[0].plot(r1,lam2)
ax[0].set_ylabel(f'Eigenvalue: $\lambda$')
ax[0].axhline(y=0,linestyle='--',color='black',linewidth=0.5)
ax[0].legend(['$\lambda_1$',f'$\lambda_2$'])

ax[1].plot(r1,W1)
ax[1].plot(r1,W2)
lim = 0.1
# ax[1].set_ylim([-lim,lim])

ax[2].plot(r1,RI[:,0,1])

plt.xlabel('Species 1')
plt.xlim([0,1]);

# Eigenvalues for reduced model are not functions of phi2
# hypothesis: outflow conditions on the boundary will be ill-posed
# for phi1 > ~0.25 because the characteristics do not move inward
# what are the characteristic variables?
