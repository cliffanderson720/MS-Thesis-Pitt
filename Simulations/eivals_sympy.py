'''
Find analytical expression for eigenvalues fo flux jacobian 
using sympy. f = [r1*(1-r1-r2)**n,
                  r2*(1-r1-r2)**n]
'''
from sympy import *
init_printing()

a1,a2,r1,r2,n,hsc,lamda = symbols('a1,a2,r1,r2,n,hsc,lambda')
hsc = (1-r1-r2)**n
flux = Matrix([a1*r1*hsc,
               a2*r2*hsc])
jac = flux.jacobian([r1,r2])
eispace_raw = jac.eigenvects()
eispace = eispace_raw.copy()

# split eigenvecs output into eigenvalues, multiplicity, and eigenvectors
# flatten the list an extract the elements
eispace = [item for component in eispace for item in component]
eivals = eispace[0::3]
multiplicity = eispace[1::3]
eivecs = eispace[2::3]
# eivecs[0] is the first eigenvector
# eivecs[0][0] is still the first eigenvector - does nothing
# eivecs[0][0][0] is the first element of the second vector
#%% Use absolute fluxes from volume balance, instead of slip velocities
flux_abs = Matrix([hsc*(a1*(1-r1)-a2*r2), # rederive this. It may not be correct
                   hsc*(a1*r1-a2*(1-r2))])
jac_abs = flux_abs.jacobian([r1,r2])
eispace_raw = jac_abs.eigenvects()
eispace = eispace_raw.copy()

eispace = [item for component in eispace for item in component]
eivals_abs = eispace[0::3]
multiplicity_abs = eispace[1::3]
eivecs_abs = eispace[2::3]


#%% Make characteristic matrix functions from sympy output
import numpy as np

def testq(a,r):
    '''Currently returns negative bacteria fluxes. Good'''
    hsc = (1-np.sum(r))**2
    cab = a*hsc-np.sum(a*hsc*r)
    qab = r*cab
    return cab,qab

def myflux_abs(a1,a2,r1,r2,n):
    '''Currently returns only positive bacterial fluxes. Not good.'''
    hsc = (1-r1-r2)**n
    return np.array([r1*hsc*(a1*(1-r1)-a2*r2),r2*hsc*(a1*r1-a2*(1-r2))])

def jacobian(a1,a2,r1,r2,n):
    jac = np.zeros((2,2))
    jac[0,0] = -a1*(-r1 - r2 + 1)**n - n*(a1*(-r1 + 1) - a2*r2)*(-r1 - r2 + 1)**n/(-r1 - r2 + 1)
    jac[0,1] =  a1*(-r1 - r2 + 1)**n - n*(a1*r1 - a2*(-r2 + 1))*(-r1 - r2 + 1)**n/(-r1 - r2 + 1)
    jac[1,0] = -a2*(-r1 - r2 + 1)**n - n*(a1*(-r1 + 1) - a2*r2)*(-r1 - r2 + 1)**n/(-r1 - r2 + 1)
    jac[1,1] =  a2*(-r1 - r2 + 1)**n - n*(a1*r1 - a2*(-r2 + 1))*(-r1 - r2 + 1)**n/(-r1 - r2 + 1)
    return jac
    
def riemann_invariants(a1,a2,r1,r2,n):
    '''Return eigenvectors of jacobian (riemann invariants)'''  
    eivecs = np.empty((2,2))
    const_term  = -(-a2*(-r1 - r2 + 1)**n - n*(a1*(-r1 + 1) - a2*r2)*(-r1 - r2 + 1)**n/(-r1 - r2 + 1))/(-a1*(-r1 - r2 + 1)**n - n*(a1*(-r1 + 1) - a2*r2)*(-r1 - r2 + 1)**n/(-r1 - r2 + 1) - (-r1 - r2 + 1)**(n - 1)*(-a1*n + a1*r1 + a1*r2 - a1 + a2*n - a2*r1 - a2*r2 + a2)/2)
    sqrt_term   = np.sqrt(-4*n*(a1**2 - 2*a1*a2 + a2**2)*(-r1 - r2 + 1)**(2*n - 1) + (-r1 - r2 + 1)**(2*n - 2)*(a1*n - a1*r1 - a1*r2 + a1 - a2*n + a2*r1 + a2*r2 - a2)**2)/2
    eivecs[0,0] = const_term + sqrt_term
    eivecs[1,0] = 1
    eivecs[0,1] = const_term - sqrt_term
    eivecs[1,1] = 1
    return eivecs

def celerity(a1,a2,r1,r2,n):
    eivals = np.empty(2)
    const_term = (-r1 - r2 + 1)**(n - 1)*(-a1*n + a1*r1 + a1*r2 - a1 + a2*n - a2*r1 - a2*r2 + a2)/2
    sqrt_term  = np.sqrt(-4*n*(a1**2 - 2*a1*a2 + a2**2)*(-r1 - r2 + 1)**(2*n - 1) + (-r1 - r2 + 1)**(2*n - 2)*(a1*n - a1*r1 - a1*r2 + a1 - a2*n + a2*r1 + a2*r2 - a2)**2)/2
    eivals[0]  = -sqrt_term + const_term
    eivals[1]  = +sqrt_term + const_term
    return eivals

states = (1,1/30,0.3,0.05,2)
myflux = myflux_abs(*states)
myjac = jacobian(*states)
ri  = riemann_invariants(*states)
ri_inv = np.linalg.inv(ri)
cel = celerity(*states)



# multiply these together to make sure they reconstruct the jacobian!!

#%% Reproducing Euler eigenvalues from the book
r,u,p,E,gamma = symbols('r,u,p,E,gamma')
E = p/(gamma-1)+1/2*r*u**2
f = Matrix([r*u,r*u**2+p,u*(E+p)])
f.jacobian([r,u*r,E])