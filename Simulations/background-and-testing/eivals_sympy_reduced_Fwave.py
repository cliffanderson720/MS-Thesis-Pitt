'''
symbolic manipulation to do eigenvalue decomposition on jacobian for a
reduced sedimentation model
'''
import sympy as sp
#%%

sp.init_printing(True)
a1,a2,a3,p1,p2,p3,n,hsc = sp.symbols('a1,a2,a3,p1,p2,p3,n,hsc')
p1l,p1r,p2l,p2r,p1_ave,p2_ave = sp.symbols('p1l,p1r,p2l,p2r,p1_ave,p2_ave')

p = sp.Matrix([p1,p2])
hsc = (1-p1)**n
flux_abs = sp.Matrix([hsc*p1*(1 -p1),
                   hsc*p2*(a2-p1)]) # for 3rd component: hsc*p3*(a3-p1)])


# getting expressions for riemann problem
pl = sp.Matrix([p1l,p2l]) # left state of riemann problem
pr = sp.Matrix([p1r,p2r]) # right state of riemann problem
# p_ave = (pl+pr)/2 # linearization of riemann problem with arithmetic average?
p_ave = sp.Matrix([p1_ave,p2_ave])
dp = pr-pl

p2pl = [(p1,p1l),(p2,p2l)]
p2pr = [(p1,p1r),(p2,p2r)]
p2p_ave = [(p1,p_ave[0]),(p2,p_ave[1])]

fluxL = flux_abs.subs(p2pl)
fluxR = flux_abs.subs(p2pr)
fluxD = fluxR-fluxL
fluxD

jac_abs = flux_abs.jacobian(p)
jac_abs

#### Characteristic decomposition
eispace_raw = jac_abs.eigenvects()
eispace = eispace_raw.copy()

# collect terms from list into eivals and eivecs
eispace = [item for component in eispace for item in component]
eivals_abs = eispace[0::3]
multiplicity_abs = eispace[1::3]
eivecs_abs = eispace[2::3]

# place list objects into sympy objects
R = sp.Matrix([[eivecs_abs[i][0][j] for i in range(len(eivecs_abs))] for j in range(len(eivecs_abs))])
Ri = R.inverse_GE()
W = Ri*p
lam = sp.Matrix(eivals_abs)

beta = Ri.subs(p2p_ave)*fluxD
beta

alpha = Ri.subs(p2p_ave)*dp
alpha


# convert sympy objects into ufuncs for evaluating with numpy
alphafunc = sp.lambdify((p1l,p1r,p1_ave,p2l,p2r,p2_ave,a2,n),alpha,'numpy')
betafunc = sp.lambdify((p1l,p1r,p1_ave,p2l,p2r,p2_ave,a2,n),beta,'numpy')
Rfunc = sp.lambdify((p1,p2,a2,n),R,'numpy')

#%% plot eigenvalues ----------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
###### Plot each as f(r2) r2 = [0.0001,0.01,0.05,0.1,0.25]
a2 = 1/30.
nr = 101
r1 = np.linspace(0,1.0,nr)
r2 = 0.1
m = 2

# eigenvalues
lam1 = ((1-r1)**m)*(a2-r1)
lam2 = ((1-r1)**m)*(-m*r1-2*r1+1)

# Riemann invariants:
RI = np.array([Rfunc(r1i,r2,a2,m) for r1i in r1])
RI[0,...]

# wave strengths computed from Ri*p
W1 = r1*(-m*r2*(a2 - r1)*(-r1 + 1)**m/(-r1 + 1) - r2*(-r1 + 1)**m)/(m*r1*(-r1 + 1)**m + r1*(-r1 + 1)**m + (a2 - r1)*(-r1 + 1)**m - (-r1 + 1)*(-r1 + 1)**m) + r2
W2 =-r1*(-m*r2*(a2 - r1)*(-r1 + 1)**m/(-r1 + 1) - r2*(-r1 + 1)**m)/(m*r1*(-r1 + 1)**m + r1*(-r1 + 1)**m + (a2 - r1)*(-r1 + 1)**m - (-r1 + 1)*(-r1 + 1)**m)

# alpha values for W waves
alphas = alphafunc(r1[:-1],r1[1:],(r1[:-1]+r1[1:])/2,r2,r2,r2,a2,m) # variable phi1, constant phi2
alpha1 = alphas[0,...].T
alpha2 = alphas[1,...].T

# beta values for F waves
betas = betafunc(r1[:-1],r1[1:],(r1[:-1]+r1[1:])/2,r2,r2,r2,a2,m)
beta1 = betas[0,...].T
beta2 = betas[1,...].T

fig,ax = plt.subplots(nrows=5,sharex=True,figsize=(10,12))
ax[0].plot(r1,lam1)
ax[0].plot(r1,lam2)
ax[0].set_ylabel('Eigenvalue: $\lambda$')
ax[0].legend(['$\lambda_1$','$\lambda_2$'])

ax[1].plot(r1,RI[:,0,0])
ax[1].plot(r1,RI[:,0,1])
lim = 1.0
# ax[1].set_ylim([-lim,lim])
ax[1].set_ylabel('$\phi_1$ in $\hat{r}^2$')
ax[1].legend(['$r^1 \phi_1$','$r^2 \phi_1$'])

ax[2].plot(r1,W1)
ax[2].plot(r1,W2)
lim = 1.0
ax[2].set_ylim([-lim,lim])
ax[2].set_ylabel('Wave strength (w)')

ax[3].plot(r1[:-1],alpha1)
ax[3].plot(r1[:-1],alpha2)
ax[3].set_ylabel('$\Delta$ wave strength ($\\alpha$)')

ax[4].plot(r1[:-1],beta1)
ax[4].plot(r1[:-1],beta2)
ax[4].set_ylabel('$\Delta$F wave strength $\\beta$')

[ax.axhline(y=0,linestyle='--',color='black',linewidth=0.5) for ax in ax]
plt.xlabel('Species 1')
plt.xlim([0,1]);

# Eigenvalues for reduced model are not functions of phi2
# hypothesis: outflow conditions on the boundary will be ill-posed
# for phi1 > ~0.25 because the characteristics do not move inward
# what are the characteristic variables?

# F waves are still discontinuous at the maximum.
