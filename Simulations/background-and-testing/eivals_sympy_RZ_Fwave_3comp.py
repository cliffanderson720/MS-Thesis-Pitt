'''
symbolic manipulation to do eigenvalue decomposition on jacobian
for the RZ zaki
'''
import sympy as sp
#%%

sp.init_printing(True)
a1,a2,a3,p1,p2,p3,n,hsc = sp.symbols('a1,a2,a3,p1,p2,p3,n,hsc')
p1l,p1r,p2l,p2r,p3l,p3r,p1_ave,p2_ave,p3_ave = sp.symbols('p1l,p1r,p2l,p2r,p3l,p3r,p1_ave,p2_ave,p3_ave')


p = sp.Matrix([p1,p2,p3])
hsc = (1-p1-p2-p3)**n
flux_abs = sp.Matrix([hsc*p1*(1 -p1-p2*a2-p3*a3),
                      hsc*p2*(a2-p1-p2*a2-p3*a3),
                      hsc*p3*(a3-p1-p2*a2-p3*a3)])

# getting expressions for riemann problem
pl = sp.Matrix([p1l,p2l,p3l]) # left state of riemann problem
pr = sp.Matrix([p1r,p2r,p3l]) # right state of riemann problem
# p_ave = (pl+pr)/2 # linearization of riemann problem with arithmetic average?
p_ave = sp.Matrix([p1_ave,p2_ave,p3_ave])
dp = pr-pl

# make dictionary-like arrays for substitution
p2pl = [(p1,p1l),(p2,p2l),(p3,p3l)]
p2pr = [(p1,p1r),(p2,p2r),(p3,p3r)]
p2p_ave = [(p1,p_ave[0]),(p2,p_ave[1]),(p3,p_ave[2])]

fluxL = flux_abs.subs(p2pl)
fluxR = flux_abs.subs(p2pr)
fluxD = fluxR-fluxL
fluxD

jac_abs = flux_abs.jacobian(p)
jac_abs

#### Characteristic decomposition
eispace_raw = jac_abs.eigenvects()
eispace = eispace_raw.copy()

# Collect terms from list into eivals and eivecs
eispace = [item for component in eispace for item in component]
eivals_abs = eispace[0::3]
multiplicity_abs = eispace[1::3]
eivecs_abs = eispace[2::3]

# convert list objects into sympy objects
R = sp.Matrix([[eivecs_abs[i][0][j] for i in range(len(eivecs_abs))] for j in range(len(eivecs_abs))])
Ri = R.inverse_LU()
W = Ri*p
lam = sp.Matrix(eivals_abs)


### Note: there is a redundancy between p and p_ave. Only 1 need exist.
beta = Ri.subs(p2p_ave)*fluxD # evaluate Ri at the interfacial value
beta

# convert sympy objects into ufuncs for evaluating with numpy
# Rfunc = sp.lambdify((p1,a2,n),R.subs(p2,0),'numpy')
# Wfunc = sp.lambdify((p1,a2,n),W.subs(p2,0),'numpy')
# lamfunc = sp.lambdify((p1,a2,n),lam.subs(p2,0),'numpy')

betafunc = sp.lambdify((p1l,p1r,p1_ave,p2l,p2r,p2_ave,a2,n),beta,'numpy')
Rfunc = sp.lambdify((p1,p2,a2,n),R,'numpy')
Rifunc = sp.lambdify((p1,p2,a2,n),Ri,'numpy')
Wfunc = sp.lambdify((p1,p2,a2,n),W,'numpy')
lamfunc = sp.lambdify((p1,p2,a2,n),lam,'numpy')
deltaFfunc = sp.lambdify((p1l,p1r,p2l,p2r,a2,n),fluxD,'numpy')


#%% plot characteristic quantities -------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

a2_= 1/30.
nr = 1001
r2 = 0.000001
r1 = np.linspace(-0.2,1-r2,nr)
dr1 = r1[1]-r1[0]
rint = (r1[:-1]+r1[1:])/2
m = 2.

# eigenvalues
lams = np.array([lamfunc(r1i,r2,a2_,m) for r1i in r1])
lam1 = lams[:,0]
lam2 = lams[:,1]

# Riemann invariants
RI = np.array([Rfunc(r1i,r2,a2_,m) for r1i in r1])
RIint = np.array([Rfunc(r1i,r2,a2_,m) for r1i in rint])

# Characteristic transform
Ris = np.array([Rifunc(r1i,r2,a2_,m) for r1i in r1])


# wave strengths computed from Ri*p
Ws = np.array([Wfunc(r1i,r2,a2_,m) for r1i in rint])
W1 = Ws[:,0]
W2 = Ws[:,1]

# beta values for F waves
betas = betafunc(r1[:-1],r1[1:],rint,r2,r2,r2,a2_,m)
beta1 = betas[0,0,:]
beta2 = betas[1,0,:]

# Style markers
species1  = dict(linestyle='-',color='blue',label='1')
species2  = dict(linestyle='-',color='red',label='2')
species12 = dict(linestyle='--',color='blue',label='12')
species21 = dict(linestyle='--',color='red',label='21')

# fig,ax = plt.subplots(nrows=6,sharex=True,figsize=(8,10))
# ax[0].plot(r1,lam1,**species1)
# ax[0].plot(r1,lam2,**species2)
# # ax[0].plot(r1,lam1+lam2)
# ax[0].set_ylabel(f'Eigenvalue: $\lambda$')
# ax[0].legend(['$\lambda_1$',f'$\lambda_2$'])
#
# ax[1].plot(r1,RI[:,0,0],**species1)
# ax[1].plot(r1,RI[:,0,1],**species12)
# ax[1].legend(['$r^1 \phi_1$','$r^2 \phi_1$'])
# lim = 2
# ax[1].set_ylim([-lim,lim])
# ax[1].set_ylabel('$\phi_1$  in R')
#
# ax[2].plot(rint,W1,**species1)
# ax[2].plot(rint,W2,**species2)
# ax[2].set_ylabel('Wave strength (w)')
# ax[2].legend(('$w_1$','$w_2$'))
#
# ax[3].plot(r1[:-1],beta1,**species1)
# ax[3].plot(r1[:-1],beta2,**species2)
# ax[3].set_ylabel('$\Delta$F wave strength $\\beta$')
#
# ax[4].plot(r1[:-1],RI[:-1,0,0]*beta1,**species1)
# ax[4].plot(r1[:-1],RI[:-1,0,1]*beta2,**species12)
# ax[4].plot(r1[:-1],RI[:-1,1,0]*beta1,**species2)
# ax[4].plot(r1[:-1],RI[:-1,1,1]*beta2,**species21)
# ax[4].set_ylabel('F waves')
#
# ax[5].plot(rint,RI[:,0,0]*W1,**species1)
# ax[5].plot(rint,RI[:,0,1]*W2,**species12)
# ax[5].plot(rint,RI[:,1,0]*W1,**species2)
# ax[5].plot(rint,RI[:,1,1]*W2,**species21)
# ax[5].set_ylabel(f'Wave strength (w)*$\\lambda$')
# ax[5].legend(('$w_1$','$w_2$'))
#
# [ax.axhline(y=0,linestyle='--',color='black',linewidth=0.5) for ax in ax]
# plt.xlabel('Species 1')
# plt.xlim([min(r1),max(r1)]);

#%%
F11 = RIint[:,0,0]*beta1
F12 = RIint[:,0,1]*beta2
F21 = RIint[:,1,0]*beta1
F22 = RIint[:,1,1]*beta2
dr1 = r1[1]-r1[0]
plt.plot(r1,lam1,**species1)
plt.plot(r1,lam2,**species12)
plt.axhline(y=0,linewidth=0.5,color='black')

fig,ax = plt.subplots(nrows=2,ncols=2,figsize=(10,10))
# plt.figure(figsize=(10,10))
ax[0,0].plot(rint,F11,**species1)
ax[0,0].plot(rint,F12,**species12)
ax[0,1].plot(rint,np.cumsum(F11)*dr1,**species1)
ax[0,1].plot(rint,np.cumsum(F12)*dr1,**species12)
ax[0,1].plot(rint,np.cumsum(F11+F12)*dr1)

ax[1,0].plot(rint,F21,**species2)
ax[1,0].plot(rint,F22,**species21)
ax[1,1].plot(rint,np.cumsum(F21)*dr1,**species2)
ax[1,1].plot(rint,np.cumsum(F22)*dr1,**species21)
ax[1,1].plot(rint,np.cumsum(F21+F22)*dr1)
[row.axhline(y=0,linewidth=0.5,linestyle='--',color='black') for a in ax for row in a]
ax[0,0].legend()
# plt.ylabel('F waves');

plt.plot(r1,r1*(1-r1)**2);plt.grid()

# RIint[:,0,1]*beta2

# Eigenvalues for reduced model are not functions of phi2
# hypothesis: outflow conditions on the boundary will be ill-posed
# for phi1 > ~0.25 because the characteristics do not move inward
# what are the characteristic variables?

# Inverse R for transformation int ocharacteristic space.
# ax[4].plot(r1,Ris[:,0,0],**species1)
# ax[4].plot(r1,Ris[:,0,1],**species12)
# ax[4].plot(r1,Ris[:,1,0],**species2)
# ax[4].plot(r1,Ris[:,1,1],**species21)
# ax[4].set_ylabel(f'Elements of $R^{-1}$')
# lim = 10
# ax[4].set_ylim([-lim,lim])
# ax[4].legend()
# ax[4].grid()
