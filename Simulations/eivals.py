'''
Find analytical expression for eigenvalues fo flux jacobian 
using sympy. f = [r1*(1-r1-r2)**n,
                  r2*(1-r1-r2)**n]
'''
from sympy import *
init_printing()

a1,a2,r1,r2,n,hsc,lamda = symbols('a1,a2,r1,r2,n,hsc,lambda')
hsc = (1-r1-r2)**n
flux = Matrix([a1*r1*hsc,a2*r2*hsc])
jac = flux.jacobian([r1,r2])
eispace_raw = jac.eigenvects()
eispace = eispace_raw.copy()

# split eigenvecs output into eigenvalues, multiplicity, and eigenvectors
# flatten the list
eispace = [item for component in eispace for item in component]
eivals = eispace[0::3]
multiplicity = eispace[1::3]
eivecs = eispace[2::3]

#%% Reproducing Euler eigenvalues from the book
r,u,p,E,gamma = symbols('r,u,p,E,gamma')
E = p/(gamma-1)+1/2*r*u**2
f = Matrix([r*u,r*u**2+p,u*(E+p)])
f.jacobian([r,u*r,E])