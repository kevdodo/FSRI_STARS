# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 15:35:12 2021

@author: kevin
"""
from math import exp
from math import sqrt

import matplotlib.pyplot as plt

day = 0
lday  = 200 

dayinc = (lday-day)/1000

L = []
T = []

class Arnett:
    def __init__(self, M, Mni):
        
        self.Mo = 2e33
        self.M = M*self.Mo
        self.Mni = Mni*self.Mo
        
        self.define_vars()
        
    def define_prim(self, M, Mni):
        self.M = M*self.Mo
        self.Mni = Mni*self.Mo
        self.define_vars()

    def define_vars(self):
        self.K = .08
        
        self.B=13.7
        self.c = 3e10
        ## 
        self.Vsc=1e9
        self.Eni0 = 4.78e10
        ## 
        
        self.Ro=9
        
        self.To = (self.K*self.M)/(self.B*self.c*self.Ro)
        self.Th = self.Ro/self.Vsc
        
        self.tm = sqrt(2*self.To*self.Th)
        self.tni =7.605e5
        
        self.eco0 = 2.561e8
        self.tco = 9.822e6
        
        self.y = (self.tm/(2*self.tni))





## co = 355 and nickel is 1##
    def D(self, x,y,s):
        tau = s*(55.3*(.1/self.K)*y**2)/((self.Vsc/1e9)*(.1+2*x*y)**2)
        g = tau/(tau+1.6)
        return g*(1+2*g*(1-g)*(1-.75*g))

    
    
    
    
    
        
    def func(self, z, y):
        return (self.D(z,y,1)*exp(-2*z*y+z**2)+5.36e-3*self.D(z,y,355)*exp(-2*.1548*z*y+z*z))*2*z
    #
    
    def simp(self,x1, y):
        
        
        x0=0
        
        n = 100
        xi = (x1-x0)/n   
        tot = 0
        st = self.func(x0,y)

        end = self.func(x1,y)
        
        
        con = 3*xi/8
        for i in range (1,n):
            x=x0+xi*i
            if i%3 ==0:
                sim =2*self.func(x, y)
            else:
                sim = 3*self.func(x,y)
            tot += sim 
        return(con*(tot+st+end))
    
    def Lambda(self, x,y):
        return exp(-x**2)*self.simp(x,y)
    
    def Lum(self, aday):
        t = 86400*aday
        xsec = t/self.tm
        #print('xsec', xsec)
        #D(xsec,y,1)*
        #return (D(xsec,y,1)*Eni0*exp(-t/tni)+
        #        D(xsec,y,355)*(eco0*exp(-t/tco)))*Mni*Lambda(xsec, y)
        return self.Lambda(xsec, self.y)*self.Mni*self.Eni0
    
        #return (Eni0*exp(-t/tni)+
        #       (eco0*exp(-t/tco)))*Mni*Lambda(xsec, y)            
    
        #return (Eni0*exp(-t/tni))*Mni*Lambda(xsec, y)

arnett = Arnett(1.4,0.6)

for i in range(8,15):
    for j in range(1,9):
        
        Mint = round(.1*i, 3)
        Mniint = round(.1*j, 3)
        arnett.define_prim(Mint, Mniint)
        
        for k in range (0 , 1001):
            day += dayinc
            L.append(arnett.Lum(day))
            
            T.append(day)
        filename="Arnett"+'M'+str(Mint) +'Mni'+ str(Mniint)
        f = open(filename+'.csv', 'w')
        f.write("time,lumin,\n")
        for t1, l1 in zip(T,L):
            output = "{:.5},{:.5},\n".format(t1, l1)
            f.write(output)
        f.close()
        L = []
        T = []
        day = 0
            
            

        
#10e43
'''

fig, ax = plt.subplots(ncols=1,nrows=1)
#ax.plot(data_x, data_y, "o")

ax.plot(T, L)

#ax.set_ylim(1e41)
ax.legend(['data1', 'data2'])
ax.set_yscale("log")

#ax.set_xscale('log')

ax.set_xlabel("time (day)")
ax.set_ylabel("luminosity (erg/s)")

'''