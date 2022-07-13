#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 08:25:55 2021

@author: emmadevin
"""

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('xtick', labelsize=14) 
matplotlib.rc('ytick', labelsize=14)
import numpy as np
import pandas as pd

x = np.linspace(-10,-4,10)
y = x




fig = plt.figure(figsize = (3.5,8))
fig.patch.set_facecolor('white')
plt.style.use('classic')


plt.subplot(311)
plt.plot(x,y, ls = '-', c = 'k', label = 'VS1')
plt.plot(x,y, ls = '--', c = 'k', label = 'VS2')
plt.plot(x,y, ls = ':', c = 'k', label = 'VS3')
plt.xscale('log')
plt.xlim(10**-3,5*10**6)
plt.ylim(1,0)
# plt.xlabel(r'Viscosity ($\eta$)', fontsize = 14)
plt.ylabel(r'depth ($D$)', fontsize = 18)
# plt.grid(which = 'both')
# plt.legend(fontsize = 12, loc='lower right')
plt.title(r'(b) VS1', fontsize =18, loc='left')
plt.xticks([10**-2,10**0,10**2,10**4,10**6])

plt.subplot(312)
plt.xscale('log')
plt.xlim(10**-3,5*10**6)
plt.ylim(1,0)
# plt.grid(which = 'both')
plt.ylabel(r'depth ($D$)', fontsize = 18)
# plt.xlabel(r'Viscosity ($\eta$)', fontsize = 14)
plt.title(r'(c) VS2', fontsize =18, loc = 'left')
plt.xticks([10**-2,10**0,10**2,10**4,10**6])

plt.subplot(313)
plt.xscale('log')
plt.xlim(10**-3,5*10**6)
# plt.grid(which = 'both')
plt.ylim(1,0)
plt.xlabel(r'viscosity ($\eta$)', fontsize = 18)
plt.ylabel(r'depth ($D$)', fontsize = 18)
plt.title(r'(d) VS3', fontsize =18, loc = 'left')
plt.xticks([10**-2,10**0,10**2,10**4,10**6])

plt.tight_layout()
plt.savefig('/Users/emmadevin/School/ViscousDissipation/Figures_FINAL/visc_grid.pdf')
plt.show()