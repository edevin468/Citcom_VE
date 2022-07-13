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


plt.subplot(231)
i = 'c01'

wav = ['', 'A','B','C','D']
box_widths = {'':2.5, 'A':1.25,'B':5.0,'C':0.625,'D':0.3125}
w = 'C'

case = '100'
FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Total_Dissipation/Case_'+case+'/CASE'+case+w+i+'.csv'

df = pd.read_csv(FILE)
l = len(df)

box_width = box_widths[w]       
t = np.array(df['Time'])
dissipation = df['Dissipation']*1


scale = 290**2 # scale for dimensionless 1m rock load
diss = dissipation/scale


plt.plot(t, diss, c ='k',ls = '-', lw = 1,label ='VS1')


case = '200'
FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Total_Dissipation/Case_'+case+'/CASE'+case+w+i+'.csv'

df = pd.read_csv(FILE)
l = len(df)

box_width = box_widths[w]       
t = np.array(df['Time'])
dissipation = df['Dissipation']*1


scale = 290**2 # scale for dimensionless 1m rock load
diss = dissipation/scale


plt.plot(t, diss, c ='k',ls = '--', lw = 2,label ='VS2')

case = '900'
FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Total_Dissipation/Case_'+case+'/CASE'+case+w+i+'.csv'

df = pd.read_csv(FILE)
l = len(df)

box_width = box_widths[w]       
t = np.array(df['Time'])
dissipation = df['Dissipation']*1


scale = 290**2 # scale for dimensionless 1m rock load
diss = dissipation/scale


plt.plot(t, diss, c ='k',ls = ':', lw = 2,label ='VS3')



plt.xlim(0,320)
plt.xlabel(r'time ($\tau_M$)', fontsize = 12)
plt.ylabel(r'energy flux ($D \mu^2 / \eta$)', fontsize = 12)
plt.title('(a)',loc ='left', fontsize = 15, y = 1.08, color = 'k')







plt.subplot(232)
i = 'c05'

wav = ['', 'A','B','C','D']
box_widths = {'':2.5, 'A':1.25,'B':5.0,'C':0.625,'D':0.3125}
w = 'C'

case = '100'
FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Total_Dissipation/Case_'+case+'/CASE'+case+w+i+'.csv'

df = pd.read_csv(FILE)
l = len(df)

box_width = box_widths[w]       
t = np.array(df['Time'])
dissipation = df['Dissipation']*1


scale = 290**2 # scale for dimensionless 1m rock load
diss = dissipation/scale


plt.plot(t, diss, c ='k',ls = '-', lw = 1,label ='VS1')


case = '200'
FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Total_Dissipation/Case_'+case+'/CASE'+case+w+i+'.csv'

df = pd.read_csv(FILE)
l = len(df)

box_width = box_widths[w]       
t = np.array(df['Time'])
dissipation = df['Dissipation']*1


scale = 290**2 # scale for dimensionless 1m rock load
diss = dissipation/scale


plt.plot(t, diss, c ='k',ls = '--', lw = 2,label ='VS2')

case = '900'
FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Total_Dissipation/Case_'+case+'/CASE'+case+w+i+'.csv'

df = pd.read_csv(FILE)
l = len(df)

box_width = box_widths[w]       
t = np.array(df['Time'])
dissipation = df['Dissipation']*1


scale = 290**2 # scale for dimensionless 1m rock load
diss = dissipation/scale


plt.plot(t, diss, c ='k',ls = ':', lw = 2,label ='VS3')


plt.xlim(0,20)
plt.xlabel(r'time ($\tau_M$)', fontsize = 12)


plt.title('(b)',loc ='left', fontsize = 15, y = 1.08, color = 'k')
leg = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), facecolor = 'w', ncol = 3, fontsize = 12)

for text in leg.get_texts():
    plt.setp(text, color = 'k')




plt.subplot(233)
i = 'c15'

wav = ['', 'A','B','C','D']
box_widths = {'':2.5, 'A':1.25,'B':5.0,'C':0.625,'D':0.3125}
w = 'C'

case = '100'
FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Total_Dissipation/Case_'+case+'/CASE'+case+w+i+'.csv'

df = pd.read_csv(FILE)
l = len(df)

box_width = box_widths[w]       
t = np.array(df['Time'])
dissipation = df['Dissipation']*1


scale = 290**2 # scale for dimensionless 1m rock load
diss = dissipation/scale


plt.plot(t, diss, c ='k',ls = '-', lw = 1,label ='VS1')


case = '200'
FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Total_Dissipation/Case_'+case+'/CASE'+case+w+i+'.csv'

df = pd.read_csv(FILE)
l = len(df)

box_width = box_widths[w]       
t = np.array(df['Time'])
dissipation = df['Dissipation']*1


scale = 290**2 # scale for dimensionless 1m rock load
diss = dissipation/scale


plt.plot(t, diss, c ='k',ls = '--', lw = 2,label ='VS2')

case = '900'
FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Total_Dissipation/Case_'+case+'/CASE'+case+w+i+'.csv'

df = pd.read_csv(FILE)
l = len(df)

box_width = box_widths[w]       
t = np.array(df['Time'])
dissipation = df['Dissipation']*1


scale = 290**2 # scale for dimensionless 1m rock load
diss = dissipation/scale


plt.plot(t, diss, c ='k',ls = ':', lw = 2,label ='VS3')

plt.xlim(0,0.02)
plt.xlabel(r'time ($\tau_M$)', fontsize = 12)
plt.title('(c)',loc ='left', fontsize = 15, y = 1.08, color = 'k')


plt.savefig('/Users/emmadevin/School/ViscousDissipation/Paper_Figures/Figure_6.eps', bbox_inches = 'tight')
# plt.tight_layout()
plt.show()
