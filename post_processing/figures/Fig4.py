#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 08:52:41 2021

@author: emmadevin
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure(figsize=(12.5, 8))
fig.patch.set_facecolor('white')
plt.style.use('classic')

case = '100'
i = 'c01'

wav = ['', 'A','B','C','D']
box_widths = {'':2.5, 'A':1.25,'B':5.0,'C':0.625,'D':0.3125} ## this is really box width not wavelength! 
w = 'C'

FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Total_Dissipation/Case_'+case+'/CASE'+case+w+i+'.csv'

df = pd.read_csv(FILE)
l = len(df)
box_width = box_widths[w]
scale = 290**2
        
t = np.array(df['Time'])
# t = [float(float_formatter(x)) for x in t]
dissipation = df['Dissipation']*1/scale
elastic = df['Elastic']*1/scale
work = df['Work']*1/scale
balance = dissipation + elastic + work


T = max(t) 


plt.subplot(231)
# plt.plot(t, balance, c ='k', ls = '-',lw = 1,label = 'Sum')
plt.plot(t, dissipation, c ='k',ls = '-', lw = 2,label ='Dissipation')
plt.plot(t, elastic, c = 'k' ,ls = '-', lw = 1, label = 'Elastic energy')
plt.plot(t, work, c='k',ls = ':', lw = 2, label = 'Work')
plt.xlim(0,320)
plt.xlabel(r'time ($\tau_M$)', fontsize = 12)
plt.ylabel(r'energy flux ($D \mu^2 / \eta$)', fontsize = 12)
plt.title('(a)',loc ='left', fontsize = 15, y = 1.08, color = 'k')



# case = '200'
i = 'c05'



FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Total_Dissipation/Case_'+case+'/CASE'+case+w+i+'.csv'

df = pd.read_csv(FILE)
l = len(df)

        
t = np.array(df['Time'])
# t = [float(float_formatter(x)) for x in t]
dissipation = df['Dissipation']*1/scale
elastic = df['Elastic']*1/scale
work = df['Work']*1/scale
balance = dissipation + elastic + work


T = max(t) 
plt.subplot(232)
# plt.plot(t, balance, c ='k', ls = '-',lw = 1,label = 'Sum')
plt.plot(t, dissipation, c ='k',ls = '-', lw = 2,label ='dissipation')
plt.plot(t, elastic, c = 'k' ,ls = '-', lw = 1, label = 'elastic energy')
plt.plot(t, work, c='k' ,ls = ':', lw = 2, label = 'work')
plt.xlim(0,20)
plt.xlabel(r'time ($\tau_M$)', fontsize = 12)


plt.title('(b)',loc ='left', fontsize = 15, y = 1.08, color = 'k')
leg = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), facecolor = 'w', ncol = 3, fontsize = 12)

for text in leg.get_texts():
    plt.setp(text, color = 'k')



# case = '500'
i = 'c15'



FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Total_Dissipation/Case_'+case+'/CASE'+case+w+i+'.csv'

df = pd.read_csv(FILE)
l = len(df)

        
t = np.array(df['Time'])
# t = [float(float_formatter(x)) for x in t]
dissipation = df['Dissipation']*1/scale
elastic = df['Elastic']*1/scale
work = df['Work']*1/scale
balance = dissipation + elastic + work


T = max(t) 
plt.subplot(233)
# plt.plot(t, balance, c ='k', ls = '-',lw = 1,label = 'Sum')
plt.plot(t, dissipation, c ='k',ls = '-', lw = 2,label ='Dissipation')
plt.plot(t, elastic, c = 'k' ,ls = '-', lw = 1, label = 'Elastic energy')
plt.plot(t, work, c='k' ,ls = ':', lw = 2, label = 'Work')
plt.xlim(0,0.02)
plt.xlabel(r'time ($\tau_M$)', fontsize = 12)
plt.title('(c)',loc ='left', fontsize = 15, y = 1.08, color = 'k')


plt.savefig('/Users/emmadevin/School/ViscousDissipation/Paper_Figures/Figure_4.eps', bbox_inches = 'tight')
# plt.tight_layout()
plt.show()
