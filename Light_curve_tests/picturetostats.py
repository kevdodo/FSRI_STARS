# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 14:26:39 2021

@author: -
"""
import pandas as pd

from statistics import median
from statistics import stdev
from statistics import mean
from math import log10, inf
import matplotlib.pyplot as plt

def analyze_data(T,L):
    max_LC = max(L)
    
    tindex = 0
    for i, l in enumerate(L):
        if l == max_LC:
            tindex = i
            break
    time = T[tindex]
    
    lday = time + 10 
    
    mindist = inf
    
    for i, t in enumerate(T):
        if t - lday < mindist and t-lday < 0:
            mindist = abs(t-lday)
            pindex = i
            
            
    print('time', T[pindex], time, lday)
    point = ((L[pindex] - L[pindex+1])/(T[pindex] - T[pindex+1])) * (lday - T[pindex]) + L[pindex]
    print('point', point, L[pindex], L[pindex+1])
    #print(L, T, pindex, point, L[pindex], L[pindex+1])
    
    try:
        
        aslope_maxL20 = abs((max_LC - point)/(time - lday))     
        
    except:
        
        aslope_maxL20 = 0.00
        
    pindex = 0
    rday = 1.5
    pointb = ((L[pindex] - L[pindex+1])/(T[pindex] - T[pindex+1])) * (rday - T[pindex]) + L[pindex]
    print(pointb, L[pindex] , L[pindex+1], max_LC, aslope_maxL20)
    
    bslope_maxL1 = (max_LC - pointb)/(time - rday)
    
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
    mad = abs(sum(MAD_LC)/N)
    
    
    
    max_LC = log10(max_LC)
    mad = log10(mad)
    aslope_maxL20 = log10(aslope_maxL20)
    bslope_maxL1 = log10(bslope_maxL1)
    median_LC = log10(median_LC)
    return max_LC, skew, CV_LC, mad, aslope_maxL20, bslope_maxL1, median_LC, time


csv_list = ["magbersten", "magstarbig", "2002ap", "2005kl", "2007gr", "sn1994I", "sn2007cl", "sn2007d", "sn2007y", "yao", "lsq"]

for i in csv_list:
    df = pd.read_csv('C:/Users/kevin/Documents/FSRI_stars/Light_curve_tests/{}.csv'.format(i))
    
    #print(df)
    T = df['time']
    L = 10**df['lumin']
    
    plt.rcParams.update({"font.size":24})
    fig, ax = plt.subplots(ncols=1,nrows=1,figsize=(12,8))
    
    ax.plot(T, L, "-o", linewidth=3, markersize=15)
    
    #ax.set_yscale("log")
    
    #ax.set_title("Light Curve: Bersten")
    ax.set_xlabel("time (day)")
    ax.set_ylabel(r"luminosity (erg s$^{-1}$)")
    ax.set_yscale("log")
    ax.set_ylim(1e41)
    
    plt.savefig(i+".jpg")
    plt.show()

    plt.close("all")

    print(i, analyze_data(T, L), '\n')
    T = []
    L = []


    
    
