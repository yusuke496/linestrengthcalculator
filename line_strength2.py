import numpy as np
from math import sqrt
from qutip import *
import numpy.linalg as LA
import matplotlib.pyplot as plt
import scipy
from scipy.integrate import quad,dblquad

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

def dm1(l0,ml0,l1,ml1):
        return complex_quadrature(lambda theta,phi:np.conjugate(scipy.special.sph_harm(ml0,l0,theta,phi))*scipy.special.sph_harm(ml1,l1,theta,phi)*np.sin(phi)*np.sin(phi)*np.exp(1j*theta)/(sqrt(2)),0,np.pi,0,2*np.pi)
def dp1(l0,ml0,l1,ml1):
        return complex_quadrature(lambda theta,phi:np.conjugate(scipy.special.sph_harm(ml0,l0,theta,phi))*scipy.special.sph_harm(ml1,l1,theta,phi)*np.sin(phi)*np.sin(phi)*np.exp(-1j*theta)/(-sqrt(2)),0,np.pi,0,2*np.pi)
def dpm0(l0,ml0,l1,ml1):
        return complex_quadrature(lambda theta,phi:np.conjugate(scipy.special.sph_harm(ml0,l0,theta,phi))*scipy.special.sph_harm(ml1,l1,theta,phi)*np.sin(phi)*np.cos(phi),0,np.pi,0,2*np.pi)
def d(j0,mj0,l0,s0,j1,mj1,l1,s1):
    dd=0+0*1j
    for ms0 in drange(-s0,s0+1,1):
        for ml0 in drange(max(-l0,mj0-ms0),min(l0,mj0-ms0)+1,1):
            for ms1 in drange(-s1,s1+1,1):
                for ml1 in drange(max(-l1,mj1-ms1),min(l1,mj1-ms1)+1,1):
                        #print(ms0,ml0,ms1,ml1)
                        #print(dp1(l0,ml0,l1,ml1)[0])
                        #print(dm1(l0,ml0,l1,ml1)[0])
                        #print(dpm0(l0,ml0,l1,ml1)[0])
                        dd+=clebsch(l0,s0,j0,ml0,ms0,mj0)*clebsch(l1,s1,j1,ml1,ms1,mj1)*dp1(l0,ml0,l1,ml1)[0]*clebsch(l1,s1,j1,ml1,ms1,mj1)*np.conjugate(dp1(l0,ml0,l1,ml1)[0])
                        dd+=clebsch(l0,s0,j0,ml0,ms0,mj0)*clebsch(l1,s1,j1,ml1,ms1,mj1)*dm1(l0,ml0,l1,ml1)[0]*clebsch(l1,s1,j1,ml1,ms1,mj1)*np.conjugate(dm1(l0,ml0,l1,ml1)[0])
                        dd+=clebsch(l0,s0,j0,ml0,ms0,mj0)*clebsch(l1,s1,j1,ml1,ms1,mj1)*dpm0(l0,ml0,l1,ml1)[0]*clebsch(l1,s1,j1,ml1,ms1,mj1)*np.conjugate(dpm0(l0,ml0,l1,ml1)[0])
    return dd
print(d(0,0,0,0,1,1,1,0))
print(dp1(0,0,1,1))
print(dm1(0,0,1,-1))
print(dpm0(0,0,1,0))
