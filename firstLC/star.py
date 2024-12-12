# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 15:35:12 2021

@author: kevin
"""
from math import exp
from math import sqrt

import matplotlib.pyplot as plt

day = 0
lday  = 100    

dayinc = (lday-day)/1000

L = []
T = []



K =.08
Mo = 2e33
M = 1.37*Mo #mej
B=13.7
c = 3e10
## 
Vsc=1e9
Eni0 = 4.78e10
## 
Mni = .6*Mo
Ro=9

To = (K*M)/(B*c*Ro)
Th = Ro/Vsc

tm = sqrt(2*To*Th)
tni =7.605e5

eco0 = 2.561e8
tco = 9.822e6

y = (tm/(2*tni))





## co = 355 and nickel is 1##
def D(x,y,s):
    tau = s*(55.3*(.1/K)*y**2)/((Vsc/1e9)*(.1+2*x*y)**2)
    g = tau/(tau+1.6)
    return g*(1+2*g*(1-g)*(1-.75*g))

def final():
    day = 1
    lday  = 100
    
    dayinc = (lday-day)/1000
    
    L = []
    T = []

    for i in range (0 , 1001):
        
        def func(z, y):

            return exp(-2*z*y+z**2)*2*z
        
        def simp(x1, y):
            
            
            x0=0
            
            n = 100
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
            return exp(-x**2)*simp(x,y)
        
        def Lum(day):
            t = 86400*day
            xsec = t/tm
            #D(xsec,y,1)*
            #return (D(xsec,y,1)*Eni0*exp(-t/tni)+
            #        D(xsec,y,355)*(eco0*exp(-t/tco)))*Mni*Lambda(xsec, y)
        
            return D(xsec,y,1)*Eni0*exp(-t/tni)*Mni*Lambda(xsec, y)+ \
                       D(xsec,y,355)*(eco0*exp(-t/tco))*Mni
        
            #return (Eni0*exp(-t/tni)+
            #       (eco0*exp(-t/tco)))*Mni*Lambda(xsec, y)            
        
            #return (Eni0*exp(-t/tni))*Mni*Lambda(xsec, y)
        
        L.append(Lum(day))
        
        T.append(day)
        day += dayinc
#10e43


    data_x = [-16.124, -15.194, -14.574, -13.953, -12.713, -12.403, -9.6124, -8.062, -5.8915, -0.31008, 3.4109, 8.3721, 14.264, 16.434, 17.674, 24.186, 29.147, 33.488, 38.76, 42.791, 47.752, 51.783, 57.674, 62.636, 67.597, 71.938, 76.899, 81.55, 87.132, 92.403, 97.364]
    data_x = [data - data_x[0] for data in data_x]
    data_y = [3.0408850256762583e+41, 4.355118736855715e+41, 7.227698036021733e+41, 1.1561122421921051e+42, 1.918668740670295e+42, 2.9991625189876286e+42, 4.87528490103389e+42, 8.016780633876855e+42, 9.817479430199784e+42, 1.2359474334445069e+43, 1.0964781961431829e+43, 8.472274141405912e+42, 5.4954087385762705e+42, 4.78630092322638e+42, 4.0926065973001263e+42, 3.2210687912834543e+42, 2.6424087573219286e+42, 2.0044720273651593e+42, 1.503141966090021e+42, 1.2217996601648811e+42, 1.0616955571987329e+42, 9.418895965228341e+41, 8.109610578538388e+41, 7.328245331389075e+41, 6.367955209079188e+41, 5.6363765582595144e+41, 5.046612975635318e+41, 4.6131757456038095e+41, 4.207266283844464e+41, 3.580964371026377e+41, 3.3036954103681355e+41]
    fig, ax = plt.subplots(ncols=1,nrows=1)
    ax.plot(data_x, data_y, "o")
    
    ax.plot(T, L)

    ax.set_ylim(1e41)
    ax.legend(['data1', 'data2'])
    ax.set_yscale("log")
    
    #ax.set_xscale('log')
    
    ax.set_xlabel("time (day)")
    ax.set_ylabel("luminosity (erg/s)")
final()

