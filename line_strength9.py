import numpy as np
from math import sqrt
from qutip import *
import numpy.linalg as LA
#import matplotlib.pyplot as plt
import scipy
import numpy
from scipy.integrate import quad,dblquad
from sympy.physics.quantum.cg import CG
import time
import tkinter

#right clic#
def make_menu(w):
    global the_menu
    the_menu = tkinter.Menu(w, tearoff=0)
    the_menu.add_command(label="cut")
    the_menu.add_command(label="copy")
    the_menu.add_command(label="paste")

def show_menu(e):
    w = e.widget
    the_menu.entryconfigure("cut", command=lambda: w.event_generate("<<Cut>>"))
    the_menu.entryconfigure("copy", command=lambda: w.event_generate("<<Copy>>"))
    the_menu.entryconfigure("paste", command=lambda: w.event_generate("<<Paste>>"))
    the_menu.tk.call("tk_popup", the_menu, e.x_root, e.y_root)
#right click

def complex_quadrature(func, a, b, c, d, **kwargs):
    def real_func(x,y):
        return numpy.real(func(x,y))
    def imag_func(x,y):
        return numpy.imag(func(x,y))
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

def ButtonEvent(event):
    ic = EditBox_nuclearspin.get()
    #ic=input("nuclear spin I:")
    i=float(ic)
    sc = EditBox_totalspin.get()
    #sc=input("total spin S:")
    s=float(sc)
    j0c = EditBox_totalJoflowerlevel.get()
    #j0c=input("total J of lower level:")
    j0=float(j0c)
    l0c = EditBox_totalloflowerlevel.get()
    #l0c=input("total L of lower level:")
    l0=float(l0c)
    j1c = EditBox_totalJofupperlevel.get()
    #j1c=input("total J of upper level:")
    j1=float(j1c)
    l1c = EditBox_totallofupperlevel.get()
    #l1c=input("total L of upper level:")
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
            
root = tkinter.Tk()
root.title("line strength")
root.geometry("400x220")

make_menu(root)
root.bind_class("Entry", "<Button-3><ButtonRelease-3>", show_menu)

Static_nuclearspin = tkinter.Label(text='nuclear spin I')
Static_nuclearspin.pack()
Static_nuclearspin.place(x=20, y=50)

Static_totalspin = tkinter.Label(text='total spin S')
Static_totalspin.pack()
Static_totalspin.place(x=20, y=70)

Static_totalJofupperlevel = tkinter.Label(text='total J of upper level')
Static_totalJofupperlevel.pack()
Static_totalJofupperlevel.place(x=20, y=90)

Static_totallofupperlevel = tkinter.Label(text='total L of upper level')
Static_totallofupperlevel.pack()
Static_totallofupperlevel.place(x=20, y=110)

Static_totalJoflowerlevel = tkinter.Label(text='total J of lower level')
Static_totalJoflowerlevel.pack()
Static_totalJoflowerlevel.place(x=20, y=130)

Static_totalloflowerlevel = tkinter.Label(text='total L of lower level')
Static_totalloflowerlevel.pack()
Static_totalloflowerlevel.place(x=20, y=150)

EditBox_nuclearspin = tkinter.Entry(width=15)
EditBox_nuclearspin.insert(tkinter.END,"4.5")
EditBox_nuclearspin.place(x=150, y=50)

EditBox_totalspin = tkinter.Entry(width=15)
EditBox_totalspin.insert(tkinter.END,"1")
EditBox_totalspin.place(x=150, y=70)

EditBox_totalJofupperlevel = tkinter.Entry(width=15)
EditBox_totalJofupperlevel.insert(tkinter.END,"3")
EditBox_totalJofupperlevel.place(x=150, y=90)

EditBox_totallofupperlevel = tkinter.Entry(width=15)
EditBox_totallofupperlevel.insert(tkinter.END,"2")
EditBox_totallofupperlevel.place(x=150, y=110)

EditBox_totalJoflowerlevel = tkinter.Entry(width=15)
EditBox_totalJoflowerlevel.insert(tkinter.END,"2")
EditBox_totalJoflowerlevel.place(x=150, y=130)

EditBox_totalloflowerlevel = tkinter.Entry(width=15)
EditBox_totalloflowerlevel.insert(tkinter.END,"1")
EditBox_totalloflowerlevel.place(x=150, y=150)
            
Button1 = tkinter.Button(text='excute', width=15)
Button1.bind("<Button-1>",ButtonEvent)#左クリック（<Button-1>）されると，ButtonEvent関数を呼び出すようにバインド
Button1.place(x=270, y=95)

Button2 = tkinter.Button(text='exit',command=root.quit, width=15)
Button2.place(x=270, y=140)

root.mainloop()
