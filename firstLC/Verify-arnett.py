# -*- coding: utf-8 -*-
"""
Created on Fri Aug 20 09:05:34 2021

@author: kevin
"""
from math import exp
from math import sqrt

day = 0
lday  = 100    

L = []
T = []



K =.08
Mo = 2e33
M = 1.4*Mo
B=13.7
c = 3e10
##
Vsc=1.2e9
Eni = 4.78e10
##
Mni = .6* Mo

Ro=9

To = (K*M)/(B*c)
Th = Ro/Vsc

tm = sqrt(2*To*Th)
tni =7.605e5

y = (tm/(2*tni))

for i in range (day , lday+1):
    
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
        return Eni*Mni *Lambda(xsec, y)
    

    L.append(Lum(day))
    
    T.append(day)
    day = day+1



print(Lambda(2, 2))

