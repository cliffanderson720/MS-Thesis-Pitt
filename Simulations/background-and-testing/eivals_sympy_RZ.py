'''
symbolic manipulation to do eigenvalue decomposition on jacobian
for the RZ zaki
'''
import sympy as sp
#%%

sp.init_printing(True)
a1,a2,a3,p1,p2,p3,n,hsc = sp.symbols('a1,a2,a3,p1,p2,p3,n,hsc')

p = sp.Matrix([p1,p2]).as_mutable()
p
hsc = (1-p1-p2)**n

flux_abs = sp.Matrix([hsc*p1*(1 -p1-p2*a2),
                      hsc*p2*(a2-p1-p2*a2)]) # for 3rd component: hsc*p3*(a3-r1)])

jac_abs = flux_abs.jacobian(p)
jac_abs

# convert to eigenspace
eispace_raw = jac_abs.eigenvects()
eispace = eispace_raw.copy()

eispace = [item for component in eispace for item in component]
eivals_abs = eispace[0::3]
multiplicity_abs = eispace[1::3]
eivecs_abs = eispace[2::3]

lam = sp.Matrix(eivals_abs)
lam
# [sp.collect(eival_abs,(-r1+1)**n) for eival_abs in eivals_abs]
# eivecs_abs

R = sp.Matrix([[eivecs_abs[i][0][j] for i in range(len(eivecs_abs))] for j in range(len(eivecs_abs))])
R
Ri = R.inverse_GE().as_mutable()

# make p2=0 for symbolic testing
W = Ri*p
W
# lim_W1 = sp.limit(W[0],p2,0)
# lim_W2 = sp.limit(W[1],p2,0)
# sp.limit(jac_abs[0,0],p2,0)
# sp.lambdify((p1,a2,n),sp.limit(lam[0],p2,0),'numpy')(0.999,1/30.,2)



# convert sympy objects into ufuncs for evaluating with numpy
# Rfunc = sp.lambdify((p1,a2,n),R.subs(p2,0),'numpy')
# Wfunc = sp.lambdify((p1,a2,n),W.subs(p2,0),'numpy')
# lamfunc = sp.lambdify((p1,a2,n),lam.subs(p2,0),'numpy')

Rfunc = sp.lambdify((p1,p2,a2,n),R,'numpy')
Wfunc = sp.lambdify((p1,p2,a2,n),W,'numpy')
lamfunc = sp.lambdify((p1,p2,a2,n),lam,'numpy')
lamfunc(r1,.0,1/30.,2)

Rfunc(r1,0.0,1/30.,2)
#%% plot characteristic quantities -------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

a2 = 1/30.
nr = 101
r2 = 0.1
r1 = np.linspace(0,1-r2,nr)
m = 3.

# eigenvalues
lams = np.array([lamfunc(r1i,r2,a2,m) for r1i in r1])
# lams = np.array([lamfunc(r1i,a2,m) for r1i in r1])
lam1 = lams[:,0]
lam2 = lams[:,1]

# Riemann invariants
RI = np.array([Rfunc(r1i,r2,a2,m) for r1i in r1])
# RI = np.array([Rfunc(r1i,a2,m) for r1i in r1])

# wave strengths computed from Ri*p
Ws = np.array([Wfunc(r1i,r2,a2,m) for r1i in r1])

# Ws = np.array([Wfunc(r1i,a2,m) for r1i in r1])
W1 = Ws[:,0]
W2 = Ws[:,1]


fig,ax = plt.subplots(nrows=3,sharex=True,figsize=(10,12))
ax[0].plot(r1,lam1)
ax[0].plot(r1,lam2)
# ax[0].plot(r1,lam1+lam2)
ax[0].set_ylabel(f'Eigenvalue: $\lambda$')
ax[0].legend(['$\lambda_1$',f'$\lambda_2$'])

ax[1].plot(r1,RI[:,0,0])
ax[1].plot(r1,RI[:,0,1])
ax[1].legend(['$r^1 \phi_1$','$r^2 \phi_1$'])
lim = 150
# ax[1].set_ylim([-lim,lim])
ax[1].set_ylabel('$\phi_1$  in R')

ax[2].plot(r1,W1)
ax[2].plot(r1,W2)
ax[2].set_ylabel('Wave strength (w)')
ax[2].legend(('$w_1$','$w_2$'))

[ax.axhline(y=0,linestyle='--',color='black',linewidth=0.5) for ax in ax]
plt.xlabel('Species 1')
plt.xlim([min(r1),max(r1)]);

# Eigenvalues for reduced model are not functions of phi2
# hypothesis: outflow conditions on the boundary will be ill-posed
# for phi1 > ~0.25 because the characteristics do not move inward
# what are the characteristic variables?



#For posterity:
#W1 = -r1/((-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m + np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2) - (-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m - np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2)) - r2*(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(((-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m + np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2) - (-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m - np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2))*(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m - np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2))

#W2 = r1/((-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m + np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2  + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2) - (-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m - np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2)) + r2*(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(((-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m + np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2) - (-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m)/(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m - np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2))*(-m*r1*(-r1 - r2 + 1)**m - r1*(-r1 - r2 + 1)**m + np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)*(-r1 - r2 + 1)**m - (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2))


# lam1 =-np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2
# lam2 = np.sqrt(4*(-r1 - r2 + 1)**(2*m)*(a2*m*r1 + a2*m*r2 + 2*a2*r1 + a2*r2 - a2 - m*r1**2 - 2*m*r1*r2 - m*r2**2 - 2*r1**2 - 4*r1*r2 + r1 - 2*r2**2 + 2*r2) + (-r1 - r2 + 1)**(2*m - 2)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)**2)/2 + (-r1 - r2 + 1)**(m - 1)*(-a2*m*r2 - a2*r1 - a2*r2 + a2 + m*r1**2 + 2*m*r1*r2 - m*r1 + m*r2**2 + 3*r1**2 + 6*r1*r2 - 4*r1 + 3*r2**2 - 4*r2 + 1)/2
