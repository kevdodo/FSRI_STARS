# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 15:37:25 2021

@author: kevin
"""

from math import exp
from math import sqrt
from math import pi
import csv
import os

import matplotlib.pyplot as plt

day = 0 #267.9
lday  = 150 #268.0
ent = 400
dayinc = (lday-day)/ent

day_loc = 0

delta = 0
n = 12
s = 0
L = []
T = []
output = {}


Bf =  1.226
Br = .987
A = .038

K =.08
Mo = 2e33
#for i in range (8,15):  
M = 0.8*Mo #mej
    
B=13.8
c = 3e10

Rosun = 6.96e10

## 

Eni0 = 4.78e10

## 
Mni = .05*Mo
Ro=9

esn = 1e51#(Vsc**2*(5/3)*M)/2

eco0 = 2.561e8
tco = 9.822e6
pcsm = 5e-13#(3*Mcsm)/(4*pi*rcsm**3)
rp = 1e14#1500*Rosun
q = pcsm*rp**s

Mcsm = .1*Mo#(4*pi/3)*rcsm*X
#change 100-3000 Ro
rcsm = (((3-s)*Mcsm/(4*pi*q))+rp**(3-s))**(1/(3-s))  #24000*Rosun
    


tocsm = (K*Mcsm/(B*c*rcsm))





#change ratio at the end 0-1
Vsc=(((10*(n-5)*esn)/(3*(n-3)*M))**(1/2))/(.5)



To = (K*M)/(B*c*Ro)
Th = Ro/Vsc
tm = sqrt(2*To*Th)
tni =7.605e5

y = (tm/(2*tni))

gn = (((1/(4*pi*(n-delta)))*(2*(5-delta)*(n-5)*esn)**((n-3)/2))/( 
    ((3-delta)*(n-3)*M)**((n-5)/2)))    



# using vsc for vsn instead saving one variable.
ti = rp/Vsc


## tfstar or trstar - t = 0 then changes function
tfstar = (((((3-s)*q**((3-n)/(n-s)))*(A*gn)**((s-3)/(n-s)))/ 
          (4*pi*Bf**(3-s)))**((n-s)/((n-3)*(3-s))))*Mcsm**((n-s)/((n-3)*(3-s)))
                
trstar = ((Vsc/(Br*(A*gn/q)**(1/(n-s))))*((1-((3-n)*M)/
                (4*pi*(Vsc**(3-n))*gn))**(1/(3-n))))**((n-s)/(s-3))


print(tfstar/86400, trstar/86400, tocsm/86400)

def thetastep(x):

    if x > 0:
        return 1
    if x<=0:
        return 0
    
    #+\
def funccsm(t):
    #print('funccsm',thetastep(tfstar-t), t)
    
    #print("{},{:.5},{:.5}".format(day,exp(t/tocsm),
    #      (t+ti)**((2*n+6*s-n*s-15)/(n-s))))
    return exp(t/tocsm)*((2*pi)/((n-s)**3))*(gn**((5-s) 
            /(n-s)))*q**((n-5)/(n-s))*((n-3)**2)*(n-5)*(Bf**(5-s))*(A**(
            (5-s)/(n-s)))*(t+ti)**((2*n+6*s-n*s-15)/(n-s))*thetastep(tfstar-t)+\
            exp(t/tocsm)*(2*pi)*((A*gn/q)**((5-n)/(n-s)))*(Br**(5-n))*gn*(((3-s)/(n-s))**3)*\
            ((t+ti)**((2*n+6*s-n*s-15)/(n-s)))*thetastep(trstar-t)
            
    
def simp(t):
            
            
    x0=0
    b = min(t, tfstar)
    n = 300
    xi = (b-x0)/n   
    tot = 0
    st = funccsm(x0)
    end = funccsm(t)
    
    
    con = 3*xi/8
    for i in range (1,n):
        x=x0+xi*i
        if i%3 ==0:
            sim =2*funccsm(x)
        else:
            sim = 3*funccsm(x)
        tot += sim 
        #print(day,x,tot)
    #print("simp", con, tot, st, end)
    return(con*(tot+st+end))
    
    
def Lambdacsm(t):
    
    return (1/tocsm)*exp(-t/tocsm)*simp(t)

#arnett
def D(x,y,s):
    tau = s*(55.3*(.1/K)*y**2)/((Vsc/1e9)*(.1+2*x*y)**2)
    g = tau/(tau+1.6)
    return g*(1+2*g*(1-g)*(1-.75*g))

def func(z, y):
    return (D(z,y,1)*exp(-2*z*y+z**2) + 5.36e-3*D(z,y,355)*exp(-2*.1548*z*y+z*z))*2*z

def simpstar(x1, y):
    
    
    x0=0
    
    n = 1000
    xi = (x1-x0)/n   
    tot = 0
    st = func(x0,y)
    end = func(x1,y)
    
    
    con = 3*xi/8
    for i in range (1,n):
        x=x0+xi*i
        if i%3 ==0:
            sim =2*func(x, y)
        else:
            sim = 3*func(x,y)
        tot += sim 
    return(con*(tot+st+end))

def Lambda(x,y):
    return exp(-x**2)*simpstar(x,y)

def Lum(day):
    t = 86400*day
    xsec = t/tm
    
    #D(xsec,y,1)*
    #return (D(xsec,y,1)*Eni0*exp(-t/tni)+
    #        D(xsec,y,355)*(eco0*exp(-t/tco)))*Mni*Lambda(xsec, y)
    return  Lambdacsm(t) #+ 
    #return Lambda(xsec, y)* Mni*Eni0
    
    #return (Eni0*exp(-t/tni)+
    #       (eco0*exp(-t/tco)))*Mni*Lambda(xsec, y)            

    #return (Eni0*exp(-t/tni))*Mni*Lambda(xsec, y)





for i in range (0 , ent):
    day_loc += dayinc
    L.append(Lum(day_loc))
    T.append(day_loc)

'''
    fig, ax = plt.subplots(ncols=1,nrows=1)
    #ax.plot(data_x, data_y, "o")
        
    ax.plot(T, L, "-")
    
    #ax.set_ylim(1e38)
    ax.set_yscale("log")
    
    #ax.set_xscale('log')
    
    ax.set_xlabel("time (day)")
    ax.set_ylabel("luminosity (erg/s)")
 
 '''
#def redefarnett(): 
    
fig, ax = plt.subplots(ncols=1,nrows=1)
#ax.plot(data_x, data_y, "o")
    
ax.plot(T, L, "-")

#ax.set_ylim(1e38)
ax.set_yscale("log")

#ax.set_xscale('log')

ax.set_xlabel("time (day)")
ax.set_ylabel("luminosity (erg/s)")

#ax.annotate(filename,(100,1e42))


    

                