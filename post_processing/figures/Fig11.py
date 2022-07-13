#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 10:45:23 2021

@author: emmadevin
"""



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure(figsize=(12.5, 6))
fig.patch.set_facecolor('white')
plt.style.use('classic')


case='900'
wav = ['', 'A','B','C','D']
D = 2.9e6
wavelengths = {'':2.5*D, 'A':1.25*D,'B':5.0*D,'C':0.625*D,'D':0.3125*D}

# data in these files in in W/m^2!!

w = 'C'

df1 = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_1/Avg_Power/Case_'+case+'/Case'+case+w+'.csv')
diss1 = df1['Dissipative_Power']
work1 = abs(df1['Work_Power'])
period1 = df1['Period']/3600/24/365



df2 = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_2/Avg_Power/Case_'+case+'/Case'+case+w+'.csv')
diss2 = df2['Dissipative_Power']
work2 = abs(df2['Work_Power'])
period2 = df2['Period']/3600/24/365



df3 = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_3/Avg_Power/Case_'+case+'/Case'+case+w+'.csv')
diss3 = df3['Dissipative_Power']
work3 = abs(df3['Work_Power'])
period3 = df3['Period']/3600/24/365




df4 = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_4/Avg_Power/Case_'+case+'/Case'+case+w+'.csv')
diss4 = df4['Dissipative_Power']
work4 = abs(df4['Work_Power'])
period4 = df4['Period']/3600/24/365


df5 = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_5/Avg_Power/Case_'+case+'/Case'+case+w+'.csv')
diss5 = df5['Dissipative_Power']
work5 = abs(df5['Work_Power'])
period5 = df5['Period']/3600/24/365


plt.subplot(121)
plt.plot(period1,diss1, c = 'k', lw=3, label=r'$\eta = 10^{22} Pa \cdot s$')
plt.plot(period1,work1, c='k', lw=3, ls = '--')

plt.plot(period2,diss2, c = 'k', lw=2.5, label=r'$\eta = 10^{21} Pa \cdot s$')
plt.plot(period2,work2, c='k', lw=2.5, ls = '--')

plt.plot(period3,diss3, c = 'k', lw=2, label=r'$\eta = 10^{20} Pa \cdot s$')
plt.plot(period3,work3, c='k', lw=2, ls = '--')

plt.plot(period4,diss4, c = 'k', lw=1.5, label=r'$\eta = 10^{19} Pa \cdot s$')
plt.plot(period4,work4, c='k', lw=1.5, ls = '--')

plt.plot(period5,diss5, c = 'k', lw=1, label=r'$\eta = 10^{18} Pa \cdot s$')
plt.plot(period5,work5, c='k', lw=1, ls = '--')

plt.text(10**-1,20**-2,'dissipation: solid line \nwork: dashed line',rotation=0)

plt.xscale('log')
plt.yscale('log')
plt.title(r'(a) $\lambda = 1.25D$',loc='left')
plt.xlabel('period (years)')
plt.ylabel(r'time-averaged energy flux ($W/m^2$)')
plt.xticks([10**-6,10**-4,10**-2,10**0,10**2,10**4,10**6])
plt.xlim(10**-4,10**6)
plt.legend(loc='upper center',bbox_to_anchor=(1.1, -0.1),ncol=5)

#########################################################
w = ''

df1 = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_1/Avg_Power/Case_'+case+'/Case'+case+w+'.csv')
diss1 = df1['Dissipative_Power']
work1 = abs(df1['Work_Power'])
period1 = df1['Period']/3600/24/365



df2 = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_2/Avg_Power/Case_'+case+'/Case'+case+w+'.csv')
diss2 = df2['Dissipative_Power']
work2 = abs(df2['Work_Power'])
period2 = df2['Period']/3600/24/365



df3 = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_3/Avg_Power/Case_'+case+'/Case'+case+w+'.csv')
diss3 = df3['Dissipative_Power']
work3 = abs(df3['Work_Power'])
period3 = df3['Period']/3600/24/365




df4 = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_4/Avg_Power/Case_'+case+'/Case'+case+w+'.csv')
diss4 = df4['Dissipative_Power']
work4 = abs(df4['Work_Power'])
period4 = df4['Period']/3600/24/365


df5 = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_5/Avg_Power/Case_'+case+'/Case'+case+w+'.csv')
diss5 = df5['Dissipative_Power']
work5 = abs(df5['Work_Power'])
period5 = df5['Period']/3600/24/365


plt.subplot(122)
plt.plot(period1,diss1, c = 'k', lw=3, label=r'$\eta = 10^{22} Pa \cdot s$')
plt.plot(period1,work1, c='k', lw=3, ls = '--')

plt.plot(period2,diss2, c = 'k', lw=2.5, label=r'$\eta = 10^{21} Pa \cdot s$')
plt.plot(period2,work2, c='k', lw=2.5, ls = '--')

plt.plot(period3,diss3, c = 'k', lw=2, label=r'$\eta = 10^{20} Pa \cdot s$')
plt.plot(period3,work3, c='k', lw=2, ls = '--')

plt.plot(period4,diss4, c = 'k', lw=1.5, label=r'$\eta = 10^{19} Pa \cdot s$')
plt.plot(period4,work4, c='k', lw=1.5, ls = '--')

plt.plot(period5,diss5, c = 'k', lw=1, label=r'$\eta = 10^{18} Pa \cdot s$')
plt.plot(period5,work5, c='k', lw=1, ls = '--')

# plt.text(10**0,20**-2,'dissipation: solid line \nwork: dashed line',rotation=0)

plt.xscale('log')
plt.yscale('log')
plt.title(r'(b) $\lambda = 5.0D$',loc='left')
plt.xlabel('period (years)')
# plt.ylabel(r'time-averaged energy flux ($W/m^2$)')
plt.xticks([10**-6,10**-4,10**-2,10**0,10**2,10**4,10**6])
plt.xlim(10**-4,10**6)





# plt.tight_layout()
# plt.legend()
plt.savefig('/Users/emmadevin/School/ViscousDissipation/Paper_Figures/Figure_11.eps', bbox_inches = 'tight')
plt.show()
