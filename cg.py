from sympy.physics.quantum.cg import CG
import numpy as np
from math import sqrt

def drange(begin, end, step):
    n = begin
    while n+step <= end:
        yield n
        n += step

def J2LS(l,s,j,jm) :
    a=[0 for i in drange(0,(2*l+1)*(2*s+1),1)]
    n=0
    for i in drange(-l,l+1,1):
        for k in drange(-s,s+1,1):
            print(k)
            a[n]=CG(l,i,s,k,j,jm).doit()
            n=n+1
            #print(a,n)
    return a

j0=1
jz0=0
s0=0.5
l0=0.5

b = J2LS(l0,s0,j0,jz0)
print(b)
