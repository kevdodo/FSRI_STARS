# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 13:35:53 2021

@author: kevin
"""
from math import exp
from math import sqrt
from math import pi
from new_star_arnett import Arnett
from statistics import median
from statistics import mean
from statistics import stdev
import matplotlib.pyplot as plt

directory = "C:/Users/kevin/Documents/FSRI_stars/csm_star/"

day = 1
lday  = 100

dayinc = (lday-day)/200

L_Arn = []
L_CSM = []
T = []
class CSM():
    def __init__(self, Mcsm, M, pcsm, rp, esn):
        self.Mo = 2e33  
        self.Rosun = 6.96e10
        
        self.Mcsm = Mcsm *self.Mo
        self.M = M*self.Mo
        self.pcsm = pcsm
        self.rp = rp
        
        self.esn = esn
        self.define_csm_vars()
        
    def define_prim(self, Mcsm, M, pcsm, rp, esn):
        #self.Mo = 2e33
        self.Mcsm = Mcsm *self.Mo
        
        self.M = M*self.Mo
        self.pcsm = pcsm
        self.rp = rp *self.Rosun
        self.esn = esn
        self.define_csm_vars()
        
    def define_csm_vars(self):
        print('rp', self.rp)
        

        self.delta = 0
        self.n = 12
        self.s = 2
        
        
        self.Bf =  1.226
        self.Br = .987
        self.A = .038
        
        self.K =.08
        #self.Mo = 2e33
        #for i in range (8,15):  
            
        self.B=13.8
        self.c = 3e10
        
        
        ## 
        
        #Eni0 = 4.78e10
    
        ## 
        #Mni = .08*Mo
        self.Ro=9
        
        #self.esn = 1e51#(Vsc**2*(5/3)*M)/2
        #self.pcsm = 5e-13
        '''
        self.pcsm = self.Mcsm/((((self.rcsm**(3-self.s))/(3-self.s))-(self.rp**(3-self.s)
                        )/(3-self.s))*4*pi*self.rp**self.s)#(3*Mcsm)/(4*pi*rcsm**3)
        '''
        self.rcsm = ((3-self.s)*((self.Mcsm/(4*pi*(self.rp**self.s)*self.pcsm)) \
                                +self.rp**(3-self.s)/(3-self.s)))**(3-self.s)
        print('rcsm' ,self.rcsm/self.Rosun)
        print('pcsm', self.pcsm)
        #self.rp = 2e14#1500*Rosun
        self.q = self.pcsm*self.rp**self.s
        
       # Mcsm = 1*self.Mo#(4*pi/3)*rcsm*X
        #change 100-3000 Ro
        #self.rcsm = (((3-self.s)*self.Mcsm/(4*pi*self.q)
        #              )+self.rp**(3-self.s))**(1/(3-self.s))  #24000*Rosun
            
        
        
        self.tocsm = (self.K*self.Mcsm/(self.B*self.c*self.rcsm))*20
        
        #(5*(self.Vsc**2)*self.M)/6
        #change ratio at the end 0-1
        self.Vsc=(((10*(self.n-5)*self.esn)/(3*(self.n-3)*self.M))**(1/2))/(.5)
        
        
        
        
        #To = (K*M)/(B*c*Ro)
        #Th = Ro/Vsc
        #tm = sqrt(2*To*Th)
        #tni =7.605e5
        
        #y = (tm/(2*tni))
        
        self.gn = (((1/(4*pi*(self.n-self.delta)))*(2*(5-self.delta)*(self.n-5)*self.esn)**((self.n-3)/2))/( 
            ((3-self.delta)*(self.n-3)*self.M)**((self.n-5)/2)))    
        
        
        
        # using vsc for vsn instead saving one variable.
        self.ti = self.rp/self.Vsc
        
        
        ## tfstar or trstar - t = 0 then changes function
        self.tfstar = (((((3-self.s)*self.q**((3-self.n)/(self.n-self.s))
                          )*(self.A*self.gn)**((self.s-3)/(self.n-self.s)))/ 
                  (4*pi*self.Bf**(3-self.s)))**((self.n-self.s)/((self.n-3)*(3-self.s)))
                       )*self.Mcsm**((self.n-self.s)/((self.n-3)*(3-self.s)))
        
        
        self.trstar = ((self.Vsc/(self.Br*(self.A*self.gn/self.q)**(1/(self.n-self.s))))*((1-((3-self.n)*self.M)/
                        (4*pi*(self.Vsc**(3-self.n))*self.gn))**(1/(3-self.n))))**((self.n-self.s)/(self.s-3))
        
        print(self.tfstar, self.trstar)
    def thetastep(self, x):
        if x > 0:
            return 1
        if x<=0:
            return 0
        #+\
    def funccsm(self, t):
    
        return exp(t/self.tocsm) * (((2*pi) / ((self.n-self.s)**3)) * (self.gn**((5-self.s) 
                /(self.n - self.s))) * self.q**((self.n - 5)/(self.n - self.s)
                ) * ((self.n - 3)**2) * (self.n - 5) * (self.Bf**(5-self.s)) * (self.A**(
                (5 - self.s) / (self.n - self.s))) * (t + self.ti)**(
                (2 * self.n + 6 * self.s - self.n * self.s-15) / (self.n - self.s)
                ) * 
                self.thetastep(self.tfstar - t) +
                (2 * pi) * ((self.A * self.gn / self.q)**(
                (5-self.n) / (self.n - self.s))) * (self.Br**(5 - self.n)) * self.gn * (
                ((3 - self.s) / (self.n - self.s))**3)*
                (t + self.ti)**((2*self.n + 6*self.s - self.n * self.s - 15)
                /(self.n - self.s)) * self.thetastep((self.trstar - t)))
        
        '''
        return ((2 * pi) * ((self.A * self.gn / self.q)**(
                (5-self.n) / (self.n - self.s))) * (self.Br**(5 - self.n)) * self.gn * (
                ((3 - self.s) / (self.n - self.s))**3))*((t + self.ti)**((2*self.n + 6*self.s - self.n * self.s - 15)
                /(self.n - self.s)) * self.thetastep((self.trstar - t)))  
                                                                         '''        


    def simp(self, t):
                
        b = min(self.trstar, t)
        x0=0
        
        n = 1000
        xi = (b-x0)/n   
        tot = 0
        st = self.funccsm(x0)
        end = self.funccsm(t)
        
        
        con = 3*xi/8
        for i in range (1,n):
            x=x0+xi*i
            if i%3 ==0:
                sim =2*self.funccsm(x)
            else:
                sim = 3*self.funccsm(x)
            tot += sim 
        return(con*(tot+st+end))
        
        
    def Lambdacsm(self, t):
        return (1/self.tocsm)*exp(-t/self.tocsm)*self.simp(t)

    def Lum(self, day):
        t = 86400*day
       # xsec = t/tm
        #D(xsec,y,1)*
        #return (D(xsec,y,1)*Eni0*exp(-t/tni)+
        #        D(xsec,y,355)*(eco0*exp(-t/tco)))*Mni*Lambda(xsec, y)
        return  self.Lambdacsm(t) #+ Lambda(xsec, y)* Mni*Eni0
        #return (Eni0*exp(-t/tni)+
        #       (eco0*exp(-t/tco)))*Mni*Lambda(xsec, y)            
    
        #return (Eni0*exp(-t/tni))*Mni*Lambda(xsec, y)
        
def analyze_data(T,L):

    max_LC = max(L)
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
    mad = sum(MAD_LC)
    
    return max_LC, skew, CV_LC, mad/N
'''
#arnett
def D(x,y,s):
    tau = s*(55.3*(.1/K)*y**2)/((Vsc/1e9)*(.1+2*x*y)**2)
    g = tau/(tau+1.6)
    return g*(1+2*g*(1-g)*(1-.75*g))

def func(z, y):
    return (D(z,y,1)*exp(-2*z*y+z**2) + 5.36e-3*D(z,y,355)*exp(-2*.1548*z*y+z*z))*2*z

def simpstar(x1, y):
    
    
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
    return exp(-x**2)*simpstar(x,y)
'''
csmstar = CSM(.1, 12, 500, 20, .5)
arnett = Arnett (12, .05)
#in solar Mass
Mcsml = [1]
Esnl = [1e52]
pcsml = [5e-6]
Rpl = [1]
# SC: Is this for Arnett or magnetar?
Ml = [4]
Mnil = [.05]
'''
f = open("LC_features.csv","w")
f.write('Mni,' + 'M,' + 'Rp,' + 'Esn,'+ 'rcsm,'+ 'mcsm,' +'max_LC,'+'skew,'+ 'CV_LC,'+ 'mad,' +'\n')
'''

for i1 in Mcsml:
    print('in iq = ', i1)
    for i2 in pcsml:
        for i3 in Esnl:
            for i4 in Rpl:
                for i5 in Ml:
                    print('in i5=', i5)
                    for i6 in Mnil: 
                        
                
                        # SC: no need to build new object, better to change the attribute 
                        # of the object
                        #magnetar = Magnetar()
                        #arnett = Arnett(1.4, .6)
                        mnili = round(i6, 3)
                        mli = round(i5, 3)
                        rpli = round(i4, 3)
                        pcsmli = round(i2, 3)
                        mcsmli = round(i1, 3)
                        esnli = round(i3, 3)
                        # SC: muted because we want to test Magnetar, not Arnett
                        arnett.define_prim(i5, i6)
                        #arnett.define_vars()
                
                        # SC: Also define define_prim for Magnetar class
                        csmstar.define_prim(i1, i5, i2, i4, i3)
                        
                        #print(i,j)
                        for i in range (0 , 201):
                            L_CSM.append(csmstar.Lum(day)) #+ 
                            L_Arn.append(arnett.Lum(day))
                            T.append(day)
                            
                            day += dayinc
                          
                        L = [l1 + l2 for l1, l2 in zip(L_CSM, L_Arn)]    
                          
                        fig, ax = plt.subplots(ncols=1,nrows=1)
                        #ax.plot(data_x, data_y, "o")
                        
                        ax.plot(T, L)
                        ax.plot(T, L_CSM)
                        
                        #ax.set_yscale("log")
                        
                        #ax.set_xscale('log')
                        
                        ax.set_xlabel("time (day)")
                        ax.set_ylabel("luminosity1 (erg/s)")
                        
                        L = []
                        L_CSM = []
                        L_Arn = []
                        '''
                        for k in range (0 , 201):
                            day += dayinc
                
                            L.append(csmstar.Lum(day)+arnett.Lum(day))
                            
                            T.append(day)
                        a, b, c, d, = analyze_data(T,L)
                        '''
                        '''
                        f.write(str(mnili) + ',' +str(mli) + ','+str(rpli) + ',' +
                      ',' + str(rcsmli) + ',' + str(mcsmli)+ ',' +
            str("{:.5}".format(a)) + ',' +str("{:.5}".format(b)) + ',' + 
            str("{:.5}".format(c)) +',' + str("{:.5}".format(d)) + '\n')
                        
                        '''
                        # SC: just for test
                        #plt.plot(T,L)
                        #plt.show()
                        
                        #filename= directory+"CSM" + "Mni" + str(mnili) + "M" + str(mli) + "rp" + str(rpli) \
                        #+ "esn" + str(esnli) + "rcsm" + str(rcsmli) + 'mcsmli' + str(rcsmli)
                        #f = open(filename+'.csv', 'w')
                        #f.write("time,lumin,\n")
                        #for t1, l1 in zip(T,L):
                        #    output = "{:.5},{:.5},\n".format(t1, l1)
                        #    f.write(output)
                        #f.close()
                        
                        #... = analyze_data(T,L)
                        
                        # format
                        # Mej,MNi,R_p,MCSM,rho_CSM,max_LC,skew,...
                                                    
       
    #f.close() 
   
                      
   
'''

for i in range (0 , 201):
    L_CSM.append(csmstar.Lum(day)) #+ 
    L_Arn.append(arnett.Lum(day))
    T.append(day)
    
    day += dayinc
  
L = [l1 + l2 for l1, l2 in zip(L_CSM, L_Arn)]    
  
fig, ax = plt.subplots(ncols=1,nrows=1)
#ax.plot(data_x, data_y, "o")

#ax.plot(T, L)
ax.plot(T, L_CSM)

#ax.set_yscale("log")

#ax.set_xscale('log')

ax.set_xlabel("time (day)")
ax.set_ylabel("luminosity1 (erg/s)")

'''