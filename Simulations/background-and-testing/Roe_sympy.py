'''
Compute Roe linearization
23 March - transfer to new file
'''

from sympy import *
init_printing(True)

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
