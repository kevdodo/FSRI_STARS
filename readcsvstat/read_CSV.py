# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 14:41:16 2021

@author: kevin
"""

import pandas as pd
from math import log10

df = pd.read_csv('C:/Users/kevin/Documents/FSRI_stars/LC_features_csm_newstat.CSV')

startype = len(df)*[2]
df['startype'] = startype

df1 = pd.read_csv('C:/Users/kevin/Documents/FSRI_stars/new_statistics/magstar_statistic_newstat.csv')

startype_mag = len(df1)*[1]
df1['startype'] = startype_mag


df2 = pd.read_csv('C:/Users/kevin/Documents/FSRI_stars/Arnettstat_newstat.csv')

startype_arn = len(df2)*[0]
df2['startype'] = startype_arn

arn_stats = df2[["max_LC", "skew", "CV_LC", "mad","aslope", "bslope", "median", "time", "startype"]]

csm_stats = df[["max_LC", "skew", "CV_LC", "mad","aslope", "bslope", "median", "time", "startype"]] 

mag_stats = df1[["max_LC", "skew", "CV_LC", "mad","aslope", "bslope", "median", "time", "startype"]]

stats = [arn_stats,  mag_stats, csm_stats]
result = pd.concat(stats)


result["mad"] = result["mad"].apply(abs).apply(log10)
result["max_LC"] = result["max_LC"].apply(log10)
result["aslope"] = result["aslope"].apply(abs)

def good_log10(num):
    if num == 0:
        return -99
    else:
        return log10(num)
    
#def add_one(x):
#    return x+1
#lambda x: x+1
    

result["aslope"] = result.apply(lambda row: good_log10(row["aslope"]), axis=1)


result["bslope"] = result["bslope"].apply(abs).apply(log10)
result["median"] = result["median"].apply(log10)
print(result.head())

f = open('C:/Users/kevin/Documents/FSRI_stars/readcsvstat/' + 'tot_lc_newstat.csv', 'w')

result.to_csv('C:/Users/kevin/Documents/FSRI_stars/readcsvstat/tot_lc_newstat.csv')
f.close()



