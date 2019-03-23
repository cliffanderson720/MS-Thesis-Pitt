'''
symbolic manipulation to find Roe parameter vector and to do
eigenvalue decomposition on jacobian
'''

from sympy import *
init_printing()

a1,a2,pl,pr,p1l,p1r,p2l,p2r,n,hsc = symbols('a1,a2,p_{l},p_{r},p_{1l},p_{1r},p_{2l},p_{2r},n,hsc')

dim = 2
wl = MatrixSymbol('wl',dim,1).as_mutable()
wr = MatrixSymbol('wr',dim,1).as_mutable()
pl = MatrixSymbol('pl',dim,1).as_mutable()
pr = MatrixSymbol('pr',dim,1).as_mutable()
B  = MatrixSymbol('B',dim,dim).as_mutable()
C  = MatrixSymbol('C',dim,dim).as_mutable()

pl[0],pl[1] = p1l,p2l
pr[0],pr[1] = p1r,p2r

wl[0],wr[0] = (sqrt(p1l),sqrt(p1r))
wl[1],wr[1] = (sqrt(p2l*p1l),sqrt(p2r*p1r))
wl,wr



def average(states,roots=False):
    '''returns symbolic mean of roots from an (n,) iterable object'''
    if roots:
        states = [sqrt(state) for state in states]
    return sum(states)/len(states)

B[0,0] = 2*average([wl[0],wr[0]])
B[0,1] = 0
B[1,0] = 0
B[1,1] = 2*average([wl[1],wr[1]])


simplify(B*(wl-wr))






#%%
a1,a2,a3,p1,p2,p3,n,hsc = symbols('a1,a2,a3,p1,p2,p3,n,hsc')

hsc = (1-p1)**n

flux_abs = Matrix([hsc*p1*(1 -p1),
                   hsc*p2*(a2-p1),
                   hsc*p3*(a3-p1)])

jac_abs = flux_abs.jacobian([p1,p2,p3])
# collect(jac_abs[1,0],(-p1+1)**n)
# collect(-p2*(1-p1)**(n-1)*(1+a2*n-n*p1-p1),(1-p1)**n)
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
[collect(eivec_abs[0][0],(-p1+1)**n) for eivec_abs in eivecs_abs]


#%%
import numpy as np
import matplotlib.pyplot as plt
a2 = 1/30.
r1 = np.linspace(0,1)
r2 = 0.00001
m = 2
lam1 = (a2-r1)*(1-r1)**m
lam2 = (1-r1)**m*(-m*r1-2*r1+1)

plt.plot(r1,lam1);plt.grid()
plt.plot(r1,lam2);plt.grid()
