import numpy as np
from math import sqrt
from qutip import *
import numpy.linalg as LA
import matplotlib.pyplot as plt
import scipy
from scipy.integrate import quad,dblquad
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j

def drange(begin, end, step):
    n = begin
    while n+step <= end:
        yield n
        n += step

def deg(fg,mfg,jg,fe,mfe,je,I):
    return (2*je+1)*(2*fg+1)*(2*fe+1)*wigner_3j(fe,1,fg,-mfe,mfe-mfg,mfg)**2*wigner_6j(fe,1,fg,jg,I,je)**2

jg=2
je=3
I=4.5
d=0
for fg in drange(I-jg,I+jg+1,1):
    for fe in drange(I-je,I+je+1,1):
        d=0
        for mfg in drange(-fg,fg+1,1):
            for mfe in drange(-fe,fe+1,1):
                #if deg(fg,mfg,jg,fe,mfe,je,I)!=0:
                #    print(mfg,"->",mfe,":",deg(fg,mfg,jg,fe,mfe,je,I))
                d+=deg(fg,mfg,jg,fe,mfe,je,I)
        if d!=0:
            print("total",fg,"->",fe,":",d)
