# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 10:34:34 2021

@author: kevin
"""

from math import exp
from math import sqrt
from math import pi
import matplotlib.pyplot as plt

from new_star_arnett import Arnett
    
day = 0
lday  = 200

dayinc = (lday-day)/1000

L = []
T = []



arnett = Arnett(1.4, .6)

#defining functions
class Magnetar():
    def __init__(self):
        self.B = 20
        self.define_vars()
    def define_vars(self):
        
        self.K =.08 # some constant
        self.Mo = 2e33 # Solar mass
        self.M = 1.4*self.Mo #mej, mass of how much gets ejected
        self.Beta=13.7 # related to density distribution 
        self.c = 3e10
        ## 
        self.Vsc=1e9
       # Eni0 = 4.78e10 # decay rate of the nickel
        ## 
        #Mni = .6*self.Mo #initial mass of nickel
       # Ro=9 # Ro is a constant, so it basically goes to zero
        
        #To = (K*self.M)/(Beta*c*Ro) 
      #  Th = Ro/Vsc 
        
       # tm = sqrt(2*To*Th)
       # tni =7.605e5
        
        #eco0 = 2.561e8 
       # tco = 9.822e6 
        
        #y = (tm/(2*tni))
        
        self.period = 0.5 # in 10 ms
        self.ep = 2e50/self.period**2
        self.tp = 1.3/(self.B**2)*86400*365*self.period**2
        
        
        
        
        self.td2 = sqrt((3/(4*pi)) * (self.M*self.K)/(self.Vsc*self.c))
        
    def func2(self, t): # Function for light curve magnetar part
    
        Lp = (1/(1+t/self.tp)**2)*(self.ep/self.tp) # magnetar luminosity - luminosity due
        #to magnetic and rotation
        
     #   print('lp' ,Lp, ep, tp)
        #t = day* 86400
       # print('func2', exp(t**2/(2*td2**2))*(Lp)*(t/td2/td2), 't', t, 'td2', td2, exp(t**2/(2*td2**2)))
        return exp(t**2/(2*self.td2**2))*(Lp)*(t/self.td2**2) 

    
    ## co = 355 and nickel is 1##

#calculates the tot rate of decay, s is the constant for ni and co


    def simp2(self, a,b): # Integral for magnetar 2, derived from differential equation
        #a and b are lower and upper bounds respectively
        n = 100
        xi = (b-a)/n   
        tot = 0
        st = self.func2(a)
        end = self.func2(b)
        
        
        con = 3*xi/8
        for i in range (1,n):
            x=a+xi*i
            if i%3 ==0:
                sim =2*self.func2(x)
            else:
                sim = 3*self.func2(x)
            tot += sim 
       # print(b, tot, con)
      #  print("IN SIMP2", con*(tot+st+end))
        return(con*(tot+st+end))

#def Lambda(x,y):
#    return exp(-x**2)*(simp(x,y)) #exp is outside lambda function
    #needs two variables for both the bounds

    def Lambda2(self,a,t): # integral and exponent part of the magnetar function
       # print("In Lambda2", b, exp(-t**2/(2*td2**2))*simp2(a,b))
    #        return #exp(-t**2/(2*td2**2))*#
        return exp(-t**2/(2*self.td2**2))*self.simp2(a,t)
    
    
    def Lum(self, day): #luminosity = magnetar + Arnett
        t = 86400*day # converts to seconds
        
        xsec = t/arnett.tm
       # return Eni0*Mni*Lambda(xsec, y)
       
        print(self.Lambda2(9,t))
        return self.Lambda2(0, t) + arnett.Eni0*arnett.Mni*arnett.Lambda(xsec, arnett.y)
    
        #return (Eni0*exp(-t/tni)+
        #       (eco0*exp(-t/tco)))*Mni*Lambda(xsec, y)            
    
        #return (Eni0*exp(-t/tni))*Mni*Lambda(xsec, y)
magnetar = Magnetar()
for i in range (0 , 1001):
#for i in range(, 101):
    
    day = i * dayinc
    L.append(magnetar.Lum(day))
    
    T.append(day)
    
#10e43

#data_x = [-16.124, -15.194, -14.574, -13.953, -12.713, -12.403, -9.6124, -8.062, -5.8915, -0.31008, 3.4109, 8.3721, 14.264, 16.434, 17.674, 24.186, 29.147, 33.488, 38.76, 42.791, 47.752, 51.783, 57.674, 62.636, 67.597, 71.938, 76.899, 81.55, 87.132, 92.403, 97.364]
#data_x = [data - data_x[0] for data in data_x]
#data_y = [3.0408850256762583e+41, 4.355118736855715e+41, 7.227698036021733e+41, 1.1561122421921051e+42, 1.918668740670295e+42, 2.9991625189876286e+42, 4.87528490103389e+42, 8.016780633876855e+42, 9.817479430199784e+42, 1.2359474334445069e+43, 1.0964781961431829e+43, 8.472274141405912e+42, 5.4954087385762705e+42, 4.78630092322638e+42, 4.0926065973001263e+42, 3.2210687912834543e+42, 2.6424087573219286e+42, 2.0044720273651593e+42, 1.503141966090021e+42, 1.2217996601648811e+42, 1.0616955571987329e+42, 9.418895965228341e+41, 8.109610578538388e+41, 7.328245331389075e+41, 6.367955209079188e+41, 5.6363765582595144e+41, 5.046612975635318e+41, 4.6131757456038095e+41, 4.207266283844464e+41, 3.580964371026377e+41, 3.3036954103681355e+41]
fig, ax = plt.subplots(ncols=1,nrows=1)
#ax.plot(data_x, data_y, "o")

ax.plot(T, L)

#ax.set_ylim(1e41)
ax.legend(['data1', 'data2'])
ax.set_yscale("log")

#ax.set_xscale('log')

ax.set_xlabel("time (day)")
ax.set_ylabel("luminosity (erg/s)")
#print (Lum(day))
#def final():
   # day = 1
   # lday  = 300
   # dayinc = (lday-day)/1000
    #t = day * 86400

'''
Extra code did not work, 

        #D(xsec,y,1)*
        #return (D(xsec,y,1)*Eni0*exp(-t/tni)+
        #        D(xsec,y,355)*(eco0*exp(-t/tco)))*Mni*Lambda(xsec, y)
      #  print("In Lum", day)
        #return  D(xsec,y,1)*Eni0*exp(-t/tni)*Mni*(Lambda(xsec, y)) +\
        #           D(xsec,y,355)*(eco0*exp(-t/tco))*Mni +\
        #               Lambda2(0, t)
'''




'''
def D(x,y,s):  #arnett
    tau = s*(55.3*(.1/K)*y**2)/((Vsc/1e9)*(.1+2*x*y)**2)
    g = tau/(tau+1.6)
    return g*(1+2*g*(1-g)*(1-.75*g))
'''

#for i in range(100, 101):


    
#def func(z,y): #function of arnett's light curve,
#D is in because it needs to be integrated because it's respect to time
#

#    return (D(z,y,1)*exp(-2*z*y+z**2) + 5.36e-3*D(z,y,355)*exp(-2*.1548*z*y+z*z))*2*z

'''
def simp(x1, y):# simpsons rule for integral arnett
    
    
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

'''
