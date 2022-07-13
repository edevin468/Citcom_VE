#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 11:16:22 2020

make work/dissipation plots

@author: emmadevin
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize = (11,8))
fig.patch.set_facecolor('white')
plt.style.use('classic')

scale = 290**2

wav = ['', 'A','B','C','D']
wavelengths = {'':2.5, 'A':1.25,'B':5.0,'C':0.625,'D':0.3125}

case='100'
df_1 = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power/Case_'+case+'/Case'+case+'C.csv')
diss_1 = df_1['Dissipative_Power']*1/scale
work_1 = abs(df_1['Work_Power'])*1/scale
period_1 = df_1['Period']

case='200'
df_2 = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power/Case_'+case+'/Case'+case+'C.csv')
diss_2 = df_2['Dissipative_Power']*1/scale
work_2 = abs(df_2['Work_Power'])*1/scale
period_2 = df_2['Period']

case='900'
df_3 = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power/Case_'+case+'/Case'+case+'C.csv')
diss_3 = df_3['Dissipative_Power']*1/scale
work_3 = abs(df_3['Work_Power'])*1/scale
period_3 = df_3['Period']



plt.subplot(221)
plt.plot(period_1,diss_1, c = 'k', lw=2.4, label='VS1')
plt.plot(period_1,work_1, c='k', lw=2.4, ls = '--')

# plt.plot(period_2,diss_2, c = 'k', ls = '-',lw=1, label='VS2')
# plt.plot(period_2,work_2, c='k', lw=1, ls = '--')

plt.plot(period_3,diss_3, c = 'k', lw=1, label='VS3')
plt.plot(period_3,work_3, c='k', lw=1, ls = '--')


    
plt.xscale('log')
plt.yscale('log')
# plt.ylim(10**-15,10**-8)

plt.xlabel(r'period ($\tau_M$)', fontsize = 12)
plt.ylabel(r'time-averaged energy flux ($ D \mu^2/\eta}$)', fontsize = 12)
plt.text(10**-2,8*10**-13,'dissipation: solid line \n work: dashed line',rotation=0)
plt.xticks([10**-5,10**-3,10**-1,10**1,10**3])
plt.title('(a)',loc='left')
plt.xlim(10**-4, 300)


plt.legend(loc='lower left')

#################
##################
#################
plt.subplot(222)
case='100'

df_C = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power/Case_'+case+'/Case'+case+'C.csv')
diss_C = df_C['Dissipative_Power']
work_C = abs(df_C['Work_Power'])
period_C = df_C['Period']


efficiency_C = diss_C/work_C

plt.plot(period_C, efficiency_C, c ='k', lw = 2.4, label='VS1')

######

case='200'

df_C = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power/Case_'+case+'/Case'+case+'C.csv')
diss_C = df_C['Dissipative_Power']
work_C = abs(df_C['Work_Power'])
period_C = df_C['Period']


efficiency_C = diss_C/work_C

# plt.plot(period_C, efficiency_C, c ='k', ls='--',lw = 1, label='VS2')

##########

case='900'


df_C = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power/Case_'+case+'/Case'+case+'C.csv')
diss_C = df_C['Dissipative_Power']
work_C = abs(df_C['Work_Power'])
period_C = df_C['Period']


efficiency_C = diss_C/work_C


plt.plot(period_C, efficiency_C, c ='k', lw = 1, label='VS3')
plt.xscale('log')
plt.yscale('log')
plt.xlim(10**-4, 300)
plt.ylim(10**-3, 2*10**0)
plt.xlabel(r'period ($\tau_M$)')
plt.ylabel('dissipation/work',labelpad=2)
plt.xticks([10**-5,10**-3,10**-1,10**1,10**3])
plt.title('(b)',loc='left')
plt.legend(loc='lower right')



plt.subplot(223)

case = '100'

df = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power/Case_'+case+'/Case'+case+'.csv')
diss = df['Dissipative_Power']*1/scale
work = abs(df['Work_Power'])*1/scale
period = df['Period']

dfA = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power/Case_'+case+'/Case'+case+'A.csv')
dissA = dfA['Dissipative_Power']*1/scale
workA = abs(dfA['Work_Power'])*1/scale
periodA = dfA['Period']

dfB = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power/Case_'+case+'/Case'+case+'B.csv')
dissB = dfB['Dissipative_Power']*1/scale
workB = abs(dfB['Work_Power'])*1/scale
periodB = dfB['Period']

dfC = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power/Case_'+case+'/Case'+case+'C.csv')
dissC = dfC['Dissipative_Power']*1/scale
workC = abs(dfC['Work_Power'])*1/scale
periodC = dfC['Period']

dfD = pd.read_csv('/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power/Case_'+case+'/Case'+case+'D.csv')
dissD = dfD['Dissipative_Power']*1/scale
workD = abs(dfD['Work_Power'])*1/scale
periodD = dfD['Period']

plt.plot(periodB,dissB, c = 'k', lw=3, label=2*wavelengths['B'])
plt.plot(periodB,workB, c='k', lw=3, ls = '--')
plt.plot(period,diss, c = 'k', lw=2, label=2*wavelengths[''])
plt.plot(period,work, c='k', lw=2, ls = '--')
# plt.plot(periodA,dissA, c = 'goldenrod', lw=1, label=2*wavelengths['A'])
# plt.plot(periodA,workA, c='goldenrod', lw=1, ls = '--')
plt.plot(periodC,dissC, c = 'k', lw=1, label=2*wavelengths['C'])
plt.plot(periodC,workC, c='k', lw=1, ls = '--')
# plt.plot(periodD,dissD, c = 'mediumblue', lw=1, label=2*wavelengths['D'])
# plt.plot(periodD,workD, c='mediumblue', lw=1, ls = '--')
plt.xscale('log')
plt.yscale('log')

plt.legend(loc='lower left', ncol = 1)
plt.title('(c)',loc='left')
plt.xlabel(r'period ($\tau_M$)', fontsize = 12)
plt.ylabel(r'time-averaged energy flux ($ D \mu^2/\eta}$)', fontsize = 12)
plt.text(10**-2,3*10**-12,'dissipation: solid line \n work: dashed line',rotation=0)
plt.xticks([10**-5,10**-3,10**-1,10**1,10**3])
plt.xlim(10**-4, 300)


plt.subplot(224)

eff = df['Dissipative_Power']/abs(df['Work_Power'])
effA = dfA['Dissipative_Power']/abs(dfA['Work_Power'])
effB = dfB['Dissipative_Power']/abs(dfB['Work_Power'])
effC = dfC['Dissipative_Power']/abs(dfC['Work_Power'])
effD = dfD['Dissipative_Power']/abs(dfD['Work_Power'])

plt.plot(period,effB, c = 'k', lw=3, label=2*wavelengths['B'])
plt.plot(period,eff, c = 'k', lw=2, label=2*wavelengths[''])
# plt.plot(period,effA, c = 'goldenrod', lw=1, label=2*wavelengths['A'])
plt.plot(period,effC, c = 'k', lw=1, label=2*wavelengths['C'])
# plt.plot(period,effD, c = 'mediumblue', lw=1, label=2*wavelengths['D'])
plt.xlim(10**-4, 300)
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='lower right', ncol = 1)
plt.title('(d)',loc='left')

plt.ylim(10**-3, 2*10**0)
plt.xlabel(r'period ($\tau_M$)')
plt.ylabel('dissipation/work',labelpad=2)
plt.xticks([10**-5,10**-3,10**-1,10**1,10**3])



plt.tight_layout()

plt.savefig('/Users/emmadevin/School/ViscousDissipation/Paper_Figures/Figure_10.eps')
plt.show()
