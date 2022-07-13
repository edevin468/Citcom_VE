#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 10:25:34 2021

@author: emmadevin
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



VELO = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE100_P/velo'
velo =  pd.read_table(VELO, sep="\s+", skiprows=1, header = None)
velo.columns = ['Node','X','Z','deltaX','deltaZ','na','na']

depth= velo['Z']
depth = depth.drop_duplicates()
depth = depth[:49]
# depth = depth[1:]
scale = 290**2

fig = plt.figure(figsize = (12,9))
plt.style.use('classic')
fig.patch.set_facecolor('white')



w = 'B'
W = 10
case = '100'

POWER = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power_by_depth/Avg_'+case+'/Case'+case+w+'.csv'
power =  pd.read_csv(POWER)
T3201 = power['320']/scale
T201 = power['20']/scale
T7631 = power['0.0195']/scale



w = 'B'
W = 10
case = '900'

POWER = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power_by_depth/Avg_'+case+'/Case'+case+w+'.csv'
power =  pd.read_csv(POWER)
T3209 = power['320']/scale
T209 = power['20']/scale
T7639 = power['0.0195']/scale

plt.subplot2grid((3,3), (0,0))

plt.plot(T3209,depth-1, c='k', lw = 1, ls = '--')
plt.plot(T3201,depth-1, c='k', lw = 1, ls = '-')

plt.ylabel(r'$y$ ($D$)', fontsize = 14) 
plt.title(r'(a) $T = 320\tau_M, \lambda = 10D$',loc ='left', fontsize = 14)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xscale('log')
plt.xlim(10**-18, 2*10**-16)

plt.subplot2grid((3,3), (1,0))

plt.plot(T201,depth-1, c='k', lw = 1, ls = '-')
plt.plot(T209,depth-1, c='k', lw = 1, ls = '--')

plt.ylabel(r'$y$ ($D$)', fontsize = 14) 
plt.title(r'(d) $T = 20\tau_M, \lambda = 10D$',loc ='left', fontsize = 14)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xscale('log')
plt.xlim(10**-16, 2*10**-14)

plt.subplot2grid((3,3), (2,0))

plt.plot(T7631,depth-1, c='k', lw = 1, ls = '-')
plt.plot(T7639, depth-1, c='k', lw = 1, ls = '--')

plt.title(r'(g) $T = 0.02\tau_M, \lambda = 10D$',loc ='left', fontsize = 14)
plt.ylabel(r'$y$ ($D$)', fontsize = 14) 
plt.xlabel('dissipation \n'+r'($\mu^2 / \eta$)', fontsize = 14)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xscale('log')
plt.xlim(10**-15, 10**-12)

w = 'A'
W = 2.5
case = '100'

POWER = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power_by_depth/Avg_'+case+'/Case'+case+w+'.csv'
power =  pd.read_csv(POWER)
T3201 = power['320']/scale
T201 = power['20']/scale
T7631 = power['0.0195']/scale

w = 'A'
W = 2.5
case = '900'

POWER = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power_by_depth/Avg_'+case+'/Case'+case+w+'.csv'
power =  pd.read_csv(POWER)
T3209 = power['320']/scale
T209 = power['20']/scale
T7639 = power['0.0195']/scale


plt.subplot2grid((3,3), (0,1))

plt.plot(T3209,depth-1, c='k', lw = 1, ls = '--')
plt.plot(T3201,depth-1, c='k', lw = 1, ls = '-')

plt.title(r'(b) $T = 320\tau_M, \lambda = 2.5D$',loc ='left', fontsize = 14)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xscale('log')
plt.xlim(10**-18, 2*10**-16)

plt.subplot2grid((3,3), (1,1))

plt.plot(T201,depth-1, c='k', lw = 1, ls = '-')
plt.plot(T209,depth-1, c='k', lw = 1, ls = '--')

plt.title(r'(e) $T = 20\tau_M, \lambda = 2.5D$',loc ='left', fontsize = 14)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xscale('log')
plt.xlim(10**-16, 2*10**-14)

plt.subplot2grid((3,3), (2,1))

plt.plot(T7631,depth-1, c='k', lw = 1, ls = '-', label='VS1')
plt.plot(T7639, depth-1, c='k', lw = 1, ls = '--', label='VS3')

plt.title(r'(h) $T = 0.02\tau_M, \lambda = 2.5D$',loc ='left', fontsize = 14)
plt.xlabel('dissipation \n'+r'($\mu^2 / \eta$)', fontsize = 14)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xscale('log')
plt.xlim(10**-15, 10**-12)
plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.7), facecolor = 'w', ncol = 2, fontsize = 14)

w = 'D'
W = 0.625
case = '100'

POWER = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power_by_depth/Avg_'+case+'/Case'+case+w+'.csv'
power =  pd.read_csv(POWER)
T3201 = power['320']/scale
T201 = power['20']/scale
T7631 = power['0.0195']/scale

w = 'D'
W = 0.625
case = '900'

POWER = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power_by_depth/Avg_'+case+'/Case'+case+w+'.csv'
power =  pd.read_csv(POWER)
T3209 = power['320']/scale
T209 = power['20']/scale
T7639 = power['0.0195']/scale

plt.subplot2grid((3,3), (0,2))

plt.plot(T3209,depth-1, c='k', lw = 1, ls = '--')
plt.plot(T3201,depth-1, c='k', lw = 1, ls = '-')

plt.title(r'(c) $T = 320\tau_M, \lambda = 0.625D$',loc ='left', fontsize = 14)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xscale('log')
plt.xticks([10**-23, 10**-21, 10**-19, 10**-17])

plt.subplot2grid((3,3), (1,2))

plt.plot(T201,depth-1, c='k', lw = 1, ls = '-')
plt.plot(T209,depth-1, c='k', lw = 1, ls = '--')

plt.title(r'(f) $T = 20\tau_M, \lambda = 0.625D$',loc ='left', fontsize = 14)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xscale('log')
plt.xticks([10**-21, 10**-18, 10**-15])

plt.subplot2grid((3,3), (2,2))

plt.plot(T7631,depth-1, c='k', lw = 1, ls = '-')
plt.plot(T7639, depth-1, c='k', lw = 1, ls = '--')

plt.title(r'(i) $T = 0.02\tau_M, \lambda = 0.625D$',loc ='left', fontsize = 14)
plt.xlabel('dissipation \n'+r'($\mu^2 / \eta$)', fontsize = 14)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xscale('log')
plt.xticks([10**-22, 10**-18, 10**-14])

plt.tight_layout()
plt.savefig('/Users/emmadevin/School/ViscousDissipation/Paper_Figures/Figure_9.eps')