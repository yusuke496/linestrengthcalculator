import numpy as np
from math import sqrt
from qutip import *
import numpy.linalg as LA
#import matplotlib.pyplot as plt
import scipy
from scipy.integrate import quad,dblquad
from sympy.physics.quantum.cg import CG
import time

def complex_quadrature(func, a, b, c, d, **kwargs):
    def real_func(x,y):
        return scipy.real(func(x,y))
    def imag_func(x,y):
        return scipy.imag(func(x,y))
    real_integral = dblquad(real_func, a, b, c, d, **kwargs)
    imag_integral = dblquad(imag_func, a, b, c, d, **kwargs)
    return (real_integral[0] + 1j*imag_integral[0], real_integral[1:], imag_integral[1:])

def drange(begin, end, step):
    n = begin
    while n+step <= end:
        yield n
        n += step

def krond(a,b):
    if a==b:
        return 1
    else:
        return 0

def d1(l0,ml0,l1,ml1,q):
        return complex_quadrature(lambda theta,phi:np.conjugate(scipy.special.sph_harm(ml0,l0,theta,phi))*scipy.special.sph_harm(ml1,l1,theta,phi)*np.sin(phi)*np.cos(phi+q*np.pi/2)*np.exp(-q*1j*theta)*np.cos(q*np.pi/4),0,np.pi,0,2*np.pi)[0]

def d2(j0,mj0,l0,s0,j1,mj1,l1,s1):
    d=0
    for ms0 in drange(-s0,s0+1,1):
        for ms1 in drange(-s1,s1+1,1):
            for ml0 in drange(max(-l0,mj0-ms0),min(l0,mj0-ms0)+1,1):
                for ml1 in drange(max(-l1,mj1-ms1),min(l1,mj1-ms1)+1,1):
                    for q in range(-1,1+1,1):
                        if ms0==ms1 and s0==s1:
                            d+=clebsch(s0,l0,j0,ms0,ml0,mj0)*clebsch(s1,l1,j1,ms1,ml1,mj1)*d1(l0,ml0,l1,ml1,q)
    return d

def d3(f0,mf0,j0,i0,l0,s0,f1,mf1,j1,i1,l1,s1):
    d=0
    for mi0 in drange(-i0,i0+1,1):
        for mi1 in drange(-i1,i1+1,1):
            for mj0 in drange(max(-j0,mf0-mi0),min(j0,mf0-mi0)+1,1):
                for mj1 in drange(max(-j1,mf1-mi1),min(j1,mf1-mi1)+1,1):
                    if mi0==mi1 and i0==i1:
                        d+=clebsch(i0,j0,f0,mi0,mj0,mf0)*clebsch(i1,j1,f1,mi1,mj1,mf1)*d2(j0,mj0,l0,s0,j1,mj1,l1,s1)
    return d

#i=4.5
#s=1
#j0=2
#l0=1
#j1=3
#l1=2

ic=input("nuclear spin I:")
i=float(ic)
sc=input("total spin S:")
s=float(sc)
j0c=input("total J of lower level:")
j0=float(j0c)
l0c=input("total L of lower level:")
l0=float(l0c)
j1c=input("total J of upper level:")
j1=float(j1c)
l1c=input("total L of upper level:")
l1=float(l1c)

for f0 in drange(np.abs(j0-i),j0+i+1,1):
    for f1 in drange(np.abs(j1-i),j1+i+1,1):
        dd=0
        end=0
        start=time.time()
        for mf0 in drange(-f0,f0+1,1) :
            for mf1 in drange(-f1,f1+1,1) :
                if np.conjugate(d3(f0,mf0,j0,i,l0,s,f1,mf1,j1,i,l1,s))*d3(f0,mf0,j0,i,l0,s,f1,mf1,j1,i,l1,s)>0.0001:
                    #print(mf0,"->",mf1,np.conjugate(d3(f0,mf0,j0,i,l0,s,f1,mf1,j1,i,l1,s))*d3(f0,mf0,j0,i,l0,s,f1,mf1,j1,i,l1,s))
                    dd+=np.conjugate(d3(f0,mf0,j0,i,l0,s,f1,mf1,j1,i,l1,s))*d3(f0,mf0,j0,i,l0,s,f1,mf1,j1,i,l1,s)
        end=time.time()-start
        if dd!=0:
            print("total",":",f0,"->",f1,":",2.5*np.real(dd),"T=",end)

input()
