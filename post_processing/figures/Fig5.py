#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 17:28:00 2021

@author: emmadevin
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

fig =plt.figure(figsize=(5,5))
fig.patch.set_facecolor('white')
plt.style.use('classic')

case = '100'
wavelengths = {'':2.5, 'A':1.25,'B':5.0,'C':0.625,'D':0.3125}
periods = [320,160,80,40,20,10,5,2.5,1.25,0.625,0.3125,0.15625,0.078125,0.0390625,0.0195,0.00976,0.00488,0.00244,0.00122,0.000610,0.000305,0.000152,0.0000763]

max_E = []
max_A = []
max_B = []
max_C = []
max_D = []

scale = 290**2


k = 1
while k < 24:
        
    if len(str(k))== 1:
        i = 'c0'+str(k)
    else:
        i = 'c'+str(k)
        
    
    k+=1

    # these are in units of D*1*mu^2/eta.
    FILE1 = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+'_P/'+i+'.time'
    FILE2 = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+'A_P/'+i+'.time'
    FILE3 = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+'B_P/'+i+'.time'
    FILE4 = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+'C_P/'+i+'.time'
    FILE5 = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+'D_P/'+i+'.time'
    
    
    df1 = pd.read_table(FILE1, sep="\s+", skiprows=1, header = None)
    df1.columns = ['Timestep','Time', 'Total Dissipation','Total Elastic','Total Elastic+Dissipation','Total Work', 'Dissipation', 'Elastic Energy','Work']
            
    df2 = pd.read_table(FILE2, sep="\s+", skiprows=1, header = None)
    df2.columns = ['Timestep','Time', 'Total Dissipation','Total Elastic','Total Elastic+Dissipation','Total Work', 'Dissipation', 'Elastic Energy','Work']
        
    df3 = pd.read_table(FILE3, sep="\s+", skiprows=1, header = None)
    df3.columns = ['Timestep','Time', 'Total Dissipation','Total Elastic','Total Elastic+Dissipation','Total Work', 'Dissipation', 'Elastic Energy','Work']
        
    df4 = pd.read_table(FILE4, sep="\s+", skiprows=1, header = None)
    df4.columns = ['Timestep','Time', 'Total Dissipation','Total Elastic','Total Elastic+Dissipation','Total Work', 'Dissipation', 'Elastic Energy','Work']
        
    df5 = pd.read_table(FILE5, sep="\s+", skiprows=1, header = None)
    df5.columns = ['Timestep','Time', 'Total Dissipation','Total Elastic','Total Elastic+Dissipation','Total Work', 'Dissipation', 'Elastic Energy','Work']
    
    
        
    t1 = np.array(df1['Time'])
    ts1 = np.array(df1['Timestep'])
    dissipation1 = list(df1['Dissipation']/wavelengths['']/scale)
    max1 = max(dissipation1)
    max_E.append(max1)
    
       
    
    t2 = np.array(df2['Time'])
    ts2 = np.array(df2['Timestep'])
    dissipation2 = list(df2['Dissipation']/wavelengths['A']/scale)
    max2 = max(dissipation2)
    max_A.append(max2)

        
    t3 = np.array(df3['Time'])
    ts3 = np.array(df3['Timestep'])
    dissipation3 = list(df3['Dissipation']/wavelengths['B']/scale)
    max3 = max(dissipation3)
    max_B.append(max3)
    
       
    t4 = np.array(df4['Time'])
    ts4 = np.array(df4['Timestep'])
    dissipation4 = list(df4['Dissipation']/wavelengths['C']/scale)
    max4 = max(dissipation4)
    max_C.append(max4)
    
      
    t5 = np.array(df5['Time'])
    ts5 = np.array(df5['Timestep'])
    dissipation5 = list(df5['Dissipation']/wavelengths['D']/scale)
    max5 = max(dissipation5)
    max_D.append(max5)


    ind1 = dissipation1.index(max1)
    ind2 = dissipation2.index(max2)
    ind3 = dissipation3.index(max3)
    ind4 = dissipation4.index(max4)
    ind5 = dissipation5.index(max5)

    # print(ts1[ind1])
    print(i, ts2[ind2])
    print(i, ts3[ind3])
    # print(ts4[ind4])
    print(i, ts5[ind5])
            
            

  
plt.plot(periods, max_B, c ='k', lw = 2, ls = '-',label=2*wavelengths['B']) 
plt.plot(periods, max_E, c ='k', lw = 2, ls = '--',label=2*wavelengths[''])
plt.plot(periods, max_A, c ='k', lw = 2, ls = '-.',label=2*wavelengths['A']) 
plt.plot(periods, max_C, c ='k', lw = 2, ls = ':',label=2*wavelengths['C']) 
plt.plot(periods, max_D, c ='k', lw = 1, ls = '-',label=2*wavelengths['D'])  

plt.xlabel(r'period ($\tau_M$)', fontsize = 12)
plt.ylabel(r'max volumetric dissipation rate ($\mu^2 / \eta$)', fontsize = 12)
plt.xscale('log')
# plt.yscale('log')
leg = plt.legend(title = r'wavelengths ($D$)',loc='upper center', bbox_to_anchor=(0.5, -0.15), facecolor = 'w', ncol = 5, fontsize = 12)

for text in leg.get_texts():
    plt.setp(text, color = 'k')


plt.savefig('/Users/emmadevin/School/ViscousDissipation/Paper_Figures/Figure_5.eps', bbox_inches='tight')
    
    
    
    
    