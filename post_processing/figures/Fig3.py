#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 10:41:37 2021

@author: emmadevin
"""


fs = 12

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

fig, (ax1,ax2,ax3) = plt.subplots(1, 3, gridspec_kw={'width_ratios': [3.3,3.3,4]},figsize=(10,4))
fig.patch.set_facecolor('white')
plt.style.use('classic')
cmap = 'gist_gray'
case = '100'
T = 'c04'
period = 40
w = 'C'
box_width = 0.625
volume = 1*1*box_width

scale = 290**2

STRESS_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'_P/'+T+'.stress.'
VELO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'_P/velo'

#====================#
# panel (a)
#====================#


ts = 60

stress = pd.read_table(STRESS_FILEPATH+str(ts), sep="\s+", skiprows=1, header = None)
stress.columns = ['SIGMA','Dissipation','ETA']
stress =stress[49:]
l = len(stress)



velo = pd.read_table(VELO_FILEPATH, sep="\s+", skiprows=1, header = None)
velo.columns = ['Node','X','Z','deltaX','deltaZ','na','na']
velo = velo[:l]  

x = velo['X']
y = velo['Z']
d = stress['Dissipation']/scale

d[-1:] = 3e-16


ax1.tricontourf(x,y-1,d, cmap = cmap,levels=20)
# cbar = plt.colorbar(format='%.e')
# cbar.set_label(r'Dissipation ($D^2 \mu^2 / \eta$)', rotation=270, labelpad=10)
ax1.set_xlabel(r'$x$ ($D$)')
ax1.set_ylabel(r'$y$ ($D$)', fontsize = 11) 
ax1.set_title(r'(a) $t = 0.8\tau_M$', loc = 'left', fontsize = fs, color = 'k')
ax1.set_xticks([0.0,0.3,0.6])


#====================#
# panel (b)
#====================#



ts = 360

stress = pd.read_table(STRESS_FILEPATH+str(ts), sep="\s+", skiprows=1, header = None)
stress.columns = ['SIGMA','Dissipation','ETA']
stress =stress[49:]
l = len(stress)



velo = pd.read_table(VELO_FILEPATH, sep="\s+", skiprows=1, header = None)
velo.columns = ['Node','X','Z','deltaX','deltaZ','na','na']
velo = velo[:l]  

x = velo['X']
y = velo['Z']
d = stress['Dissipation']/scale

d[-1:] = 3e-16


ax2.tricontourf(x,y-1,d , cmap = cmap, levels=20)
# cbar = plt.colorbar(format='%.e')
# cbar.set_label(r'Dissipation ($D^2 \mu^2 / \eta$)', rotation=270, labelpad=10)
ax2.set_xlabel(r'$x$ ($D$)')
ax2.set_yticklabels([])
# plt.ylabel(r'y ($D$)', fontsize = 11) 
ax2.set_title(r'(b) $t = 5\tau_M$', loc = 'left', fontsize = fs, color = 'k')
ax2.set_xticks([0.0,0.3,0.6])



#====================#
# panel (c)
#====================#


ts = 780

stress = pd.read_table(STRESS_FILEPATH+str(ts), sep="\s+", skiprows=1, header = None)
stress.columns = ['SIGMA','Dissipation','ETA']
stress =stress[49:]
l = len(stress)



velo = pd.read_table(VELO_FILEPATH, sep="\s+", skiprows=1, header = None)
velo.columns = ['Node','X','Z','deltaX','deltaZ','na','na']
velo = velo[:l]  

x = velo['X']
y = velo['Z']
d = stress['Dissipation']/scale

d[-1:] = 3e-16


im = ax3.tricontourf(x,y-1,d, cmap = cmap, levels=20)
cbar = fig.colorbar(im,format='%.1e')
cbar.set_label(r'dissipation ( $\mu^2 / \eta$)', rotation=270, labelpad=20)
ax3.set_xlabel(r'$x$ ($D$)')
ax3.set_yticklabels([])
# plt.ylabel(r'y ($D$)', fontsize = 11) 
ax3.set_title(r'(c) $t = 10\tau_M$', loc = 'left', fontsize = fs, color = 'k')
ax3.set_xticks([0.0,0.3,0.6])


   


#====================#
# color bar
#====================#


# plt.subplot(144)
# cbar = plt.colorbar(format='%.e')
# cbar.set_label(r'Dissipation ($D^2 \mu^2 / \eta$)', rotation=270, labelpad=20)


plt.tight_layout()
plt.savefig('/Users/emmadevin/School/ViscousDissipation/Paper_Figures/Figure_3.eps', bbox_inches = 'tight')

