# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 10:34:34 2021

@author: kevin
"""
import sys
from math import exp
from math import sqrt
from math import pi

from statistics import median
from statistics import mean
from statistics import stdev



sys.path.append('C:/Users/kevin/Documents/FSRI_stars/')

from new_star_arnett import Arnett
import matplotlib.pyplot as plt


from statistics import quantiles

import numpy as np

day = 1
lday  = 100

dayinc = (lday-day)/200

L = []
T = []
L_mag = []
L_Arn = []
#defining functions
class Magnetar():

    # SC: Better pass B as an input parameter
    def __init__(self, M, B, P, Mniint):
        self.Mo = 2e33
        self.B = B
        self.period = P 
        self.M = M * self.Mo
        self.Mni = Mniint *self.Mo
        self.define_vars()
        
        
    def define_prim(self, M, B, P, Mniint):
        self.Mo = 2e33
        self.B = B
        self.period = P 
        self.M = M *self.Mo
        
        self.Mni = Mniint * self.Mo
        self.define_vars()
        
    def define_vars(self):
        
        self.K =.08 # some constant
        self.Mo = 2e33 # Solar mass

        # SC: this is also a model parameter so better set it as parameter
        # to be passed when this method is called 
        #mej, mass of how much gets ejected

        self.Beta=13.7 # related to density distribution 
        self.c = 3e10
        ## 
        self.Vsc= sqrt(2.53*(self.Mni/self.M))*1e9 #1e9
        self.Esn = (5*(self.Vsc**2)*self.M)/6
        
        
        #Eni0 = 4.78e10 # decay rate of the nickel
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
        
        
        
        # SC: same here that it should be passed as a parameter,
        # not set by hand 
        # p in 10 ms

        self.ep = 2e50/self.period**2
        self.tp = 1.3/(self.B**2)*86400*365*self.period**2
        self.vf = sqrt((self.ep + self.Esn)/(2*self.M))
        self.td2 = sqrt((3/(4*pi)) * (self.M*self.K)/(self.vf*self.c))
        
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
        n = 500
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

    # SC: The notation a, t is consfusing 
    def Lambda2(self,a,t): # integral and exponent part of the magnetar function
       # print("In Lambda2", b, exp(-t**2/(2*td2**2))*simp2(a,b))
    #        return #exp(-t**2/(2*td2**2))*#
        return exp(-t**2/(2*self.td2**2))*self.simp2(a,t) 
    
    
    def Lum(self, day): #luminosity = magnetar + Arnett
        t = 86400*day # converts to seconds
       # return Eni0*Mni*Lambda(xsec, y)
       
       
        # SC: we don't need to output arnett inside magnetar
        return self.Lambda2(0, t)# + arnett.Lum(t)
    
        #return (Eni0*exp(-t/tni)+
        #       (eco0*exp(-t/tco)))*Mni*Lambda(xsec, y)            

        #return (Eni0*exp(-t/tni))*Mni*Lambda(xsec, y)
        
def analyze_data(T,L):
    
    
    #print(first_quartile)
    max_LC = max(L)
    
    tindex = 0
    for i, l in enumerate(L):
        if l == max_LC:
            tindex = i
            break
    time = T[tindex]
    
    try:
        
        aslope_maxL20 = (max_LC - L[20+tindex])/(T[tindex] - (10 + time))
        
    except:
        
        aslope_maxL20 = 0.00
        
    bslope_maxL1 = (max_LC - L[1])/(T[tindex] - .5)
    
    median_LC = median(L)
    
    skewl = []
    MAD_LC = []
    
    N = len(L)
    
    for i in L:
        j = (i - mean(L))**4/N
        skewl.append(j)
    skew = sum(skewl)/stdev(L)**4
    CV_LC = stdev(L)/mean(L)
    
    for k in L:
        MAD = k - median_LC
        MAD_LC.append(MAD)
    mad = sum(MAD_LC)/N
    
    return max_LC, skew, CV_LC, mad, aslope_maxL20, bslope_maxL1, median_LC, time

magnetar = Magnetar(1.4, .5, .2, .5)
arnett = Arnett(1.4, .6)



f = open('magstar_statistic_newstat.CSV', 'w')
f.write('B,' + 'p10,' + 'Mint,' + 'Mniint,'+ 'max_LC,'+'skew,'+ 'CV_LC,'+ 'mad,' 
        + 'aslope,' + 'bslope,' + 'median,' + 'time' +'\n')



p10 = [.2, .4, .6, .8, 1, 1.1, 1.4, 1.6, 1.8, 2, 10, 12, 15, 17, 20, 200]
    
B = [.5, 1, 2, 3, 5.5, 10, 12, 15, 20, 22.5, 25, 30, 35]

mass_list = [.4, .5, 1, 2, 3, 5, 7, 10]

mni_list = [.1, 0.3, .4, 0.5, 0.6, 0.8]

#print(len(p10), len(B), len(mass_list), len(mni_list))
# SC: Is this for Arnett or magnetar?
for q in p10: #p10:
    print("In q = ", q)
    for a1 in B:
        for i in mass_list:
            for j in mni_list:
        
                # SC: no need to build new object, better to change the attribute 
                # of the object
                #magnetar = Magnetar()
                #arnett = Arnett(1.4, .6)
                qint = q
                aint = a1
                Mint = i
                Mniint = j
                
                # SC: muted because we want to test Magnetar, not Arnett
                arnett.define_prim(Mint, Mniint)
                #arnett.define_vars()
        
                # SC: Also define define_prim for Magnetar class
                magnetar.define_prim(Mint, a1, q, Mniint)
                
                #print(i,j)
                for k in range (0 , 201):
                    day += dayinc
        
                    L_mag.append(magnetar.Lum(day))
                    L_Arn.append(arnett.Lum(day))
                    
                    T.append(day)
                L = [l1 + l2 for l1, l2 in zip(L_mag, L_Arn)] 
                
                '''
                fig, ax = plt.subplots(ncols=1,nrows=1)
                #ax.plot(data_x, data_y, "o")
                
                #ax.plot(T, L)
            
                ax.plot(T, L_mag)
                
                #ax.set_yscale("log")
                
                #ax.set_xscale('log')
                
                ax.set_xlabel("time (day)")
                ax.set_ylabel("luminosity1 (erg/s)")
                '''
                
                A, B1, C, D, E, F, G , H = analyze_data(T,L)
                
                #f.write('p10 ' + str(qint) + 'B ' + str(aint) + "M " + str(Mint) +
                #             "Mni " + str(Mniint) )
                
                f.write(str(aint) + ',' + str(qint) + ',' + str(Mint)+ ',' +
                str(Mniint) + ',' + str("{:.5}".format(A)) + ',' +str("{:.5}".format(B1)) + ',' + 
                str("{:.5}".format(C)) +',' + str("{:.5}".format(D)) + ',' + 
                str("{:.5}".format(E)) + ',' + str("{:.5}".format(F)) + ',' + 
                str("{:.5}".format(G)) + ',' + str("{:.5}".format(H)) +'\n')
                
                L_mag = []
                L_Arn = []
                L = []
                T = []    
                day = 0
                # SC: just for test
                #plt.plot(T,L)
                #plt.show()
                
                
                
                
'''
filename="Magnetar"+'M'+str(Mint) +'Mni'+ str(Mniint) + 'qint' + str(qint) + 'aint' + str(aint)
f = open('C:/Users/kevin/Documents/FSRI_stars/magnater_csv' + filename+'.csv', 'w')
f.write("time,lumin,\n")
for t1, l1 in zip(T,L):
output = "{:.5},{:.5},\n".format(t1, l1)
f.write(output)
'''