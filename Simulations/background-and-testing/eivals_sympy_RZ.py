'''
symbolic manipulation to do eigenvalue decomposition on jacobian
for the RZ zaki
'''
from sympy import *
#%%

init_printing(False)
a1,a2,a3,r1,p2,p3,n,hsc = symbols('a1,a2,a3,r1,p2,p3,n,hsc')
p = Matrix([r1,p2]).as_mutable()
p
hsc = (1-r1-p2)**n

flux_abs = Matrix([hsc*r1*(1 -r1-p2),
                   hsc*p2*(a2-r1-p2)]) # for 3rd component: hsc*p3*(a3-r1)])

jac_abs = flux_abs.jacobian(p)
jac_abs

# convert to eigenspace
eispace_raw = jac_abs.eigenvects()
eispace = eispace_raw.copy()

eispace = [item for component in eispace for item in component]
eivals_abs = eispace[0::3]
multiplicity_abs = eispace[1::3]
eivecs_abs = eispace[2::3]

[collect(eival_abs,(-r1+1)**n) for eival_abs in eivals_abs]
eivecs_abs

R = Matrix([[eivecs_abs[i][0][j] for i in range(len(eivecs_abs))] for j in range(len(eivecs_abs))]).as_mutable()
Ri = R.inverse_GE().as_mutable()
W=Ri*p

Rfunc = lambdify((p1,p2,a2,n),R,'numpy')


#%% plot eigenvalues
import numpy as np
import matplotlib.pyplot as plt
a2 = 1/30.
nr = 101
r1 = np.linspace(0.01,0.95,nr)
r2 = 0.01
m = 2

# eigenvalues
lam1 =-np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2
lam2 = np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2


# waves computed from Ri*p
W1 = -r1/((-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m + np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2) - (-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m - np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2)) - r2*(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(((-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m + np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2) - (-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m - np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2))*(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m - np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2))

W2 = r1/((-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m + np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2  + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2) - (-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m - np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2)) + r2*(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(((-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m + np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2) - (-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m - np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2))*(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m + np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2))

# Riemann invariants
RI = np.array([Rfunc(r1i,r2,a2,m) for r1i in r1])
Rfunc(0.2,0.01,1/30.,2)

fig,ax = plt.subplots(nrows=2,sharex=True)
ax[0].plot(r1,lam1)
ax[0].plot(r1,lam2)
ax[0].set_ylabel(f'Eigenvalue: $\lambda$')
ax[0].axhline(y=0,linestyle='--',color='black',linewidth=0.5)
ax[0].legend(['$\lambda_1$',f'$\lambda_2$'])

ax[1].plot(r1,RI2)
# ax[1].plot(r1,)
lim = .1
# ax[1].set_ylim([-lim,lim])

plt.xlabel('Species 1')
plt.xlim([0,1]);

# Eigenvalues for reduced model are not functions of phi2
# hypothesis: outflow conditions on the boundary will be ill-posed
# for phi1 > ~0.25 because the characteristics do not move inward
# what are the characteristic variables?
