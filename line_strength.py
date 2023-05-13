#import numpy as np
#from math import sqrt
#from qutip import *
#import numpy.linalg as LA
#import matplotlib.pyplot as plt
#import scipy

def drange(begin, end, step):
    n = begin
    while n+step <= end:
        yield n
        n += step
def Fm1(f,j,i):
    return (f+j-i-1)*(f+j-i)*(i+j+f+1)*(i+j+f)/4/f

def F0(f,j,i):
    return (2*f+1)*(f+j-i)*(f+i-j+1)*(i+j+1+f)*(i+j-f)/4/f/(f+1)

def Fp1(f,j,i):
    return (f+i-j+1)*(j+i-f)*(f+i-j+2)*(j+i-f-1)/4/(f+1)

jj=float(input("j="))
ii=float(input("i="))
for ff in drange(abs(jj-ii),jj+ii+1,1):
    for i0 in range(-1,2):
        if i0==-1:
            if Fm1(ff,jj,ii)!=0:
                print(ff,"->",ff-1,":",Fm1(ff,jj,ii))
        elif i0==0:
            if F0(ff,jj,ii)!=0:
                print(ff,"->",ff,":",F0(ff,jj,ii))
        elif i0==1:
            if Fp1(ff,jj,ii)!=0:
                print(ff,"->",ff+1,":",Fp1(ff,jj,ii))
