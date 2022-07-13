#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 10:41:37 2021

THIS ONE IS RIGHT

@author: emmadevin
"""




import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

fig = plt.figure(figsize = (12,9))

fig.patch.set_facecolor('white')
plt.style.use('classic')
cmap = 'gist_gray'

case = '900'

scale = 290**2

T = 'c01'
period = 320
wav = ['', 'A','B','C','D']
fs = 13

w = wav[2]
TOPO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.topo_s_'
STRESS_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.stress.'
VELO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/velo'

TIME_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.time'
time = pd.read_table(TIME_FILEPATH, sep="\s+", skiprows=1, header = None)
time.columns = ['Timestep','Time', 'Total Dissipation','Total Elastic',
                    'Total Elastic+Dissipation','Total Work', 'Dissipation', 
                    'Elastic Energy','Work']

diss = time['Dissipation']
max_diss = max(diss)
ts = time.loc[time['Dissipation'] == max_diss, 'Timestep'].iloc[0] 

if len(str(ts)) == 3:  t = '00'+str(ts)
elif len(str(ts)) == 2:  t = '000'+str(ts)
elif len(str(ts)) == 1:  t = '0000'+str(ts)



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

# d[-1:] = 4.5e-10


    
plt.subplot2grid((3,3), (0,0))  
plt.tricontourf(x,y-1,d, cmap=cmap,levels=20)
# cbar = plt.colorbar(format='%.1e')
# cbar.set_label(r'Dissipation ($D^2 \mu^2 / \eta$)', rotation=270, labelpad=10)
# plt.xlabel(r'Horizontal   $x$ ($D$)', fontsize = 14) 
plt.ylabel(r'  $y$ ($D$)', fontsize = 14) 
plt.xticks([])
plt.title(r'(a) $T = 320\tau_M$, $\lambda = 10D$',loc ='left', fontsize = fs)
cbar = plt.colorbar(format='%.1e')
 
  #####



  #####
w = wav[1]
TOPO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.topo_s_'
STRESS_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.stress.'
VELO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/velo'

TIME_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.time'
time = pd.read_table(TIME_FILEPATH, sep="\s+", skiprows=1, header = None)
time.columns = ['Timestep','Time', 'Total Dissipation','Total Elastic',
                    'Total Elastic+Dissipation','Total Work', 'Dissipation', 
                    'Elastic Energy','Work']

diss = time['Dissipation']
max_diss = max(diss)
ts = time.loc[time['Dissipation'] == max_diss, 'Timestep'].iloc[0]+150

if len(str(ts)) == 5:  t = str(ts)
elif len(str(ts)) == 4:  t = '0'+str(ts)
elif len(str(ts)) == 3:  t = '00'+str(ts)
elif len(str(ts)) == 2:  t = '000'+str(ts)
elif len(str(ts)) == 1:  t = '0000'+str(ts)



# read in relavent files for each time step i
topo = pd.read_table(TOPO_FILEPATH+str(t)+'.dat', sep="\s+", skiprows=1, header = None)
topo.columns = ['x','z','zprime','na','na']

stress = pd.read_table(STRESS_FILEPATH+str(ts), sep="\s+", skiprows=1, header = None)
stress.columns = ['SIGMA','Dissipation','ETA']
stress =stress[49:]
l = len(stress)

# visc = pd.read_table(VISC_FILEPATH, sep="\s+", skiprows=1, header = None)
# visc.columns = ['na','na','visc']

velo = pd.read_table(VELO_FILEPATH, sep="\s+", skiprows=1, header = None)
velo.columns = ['Node','X','Z','deltaX','deltaZ','na','na']
velo = velo[:l]  


x = velo['X']
y = velo['Z']
d = stress['Dissipation']/scale
# d[-1:] = 4.5e-10


    
plt.subplot2grid((3,3), (0,1))  
plt.tricontourf(x,y-1,d, cmap=cmap,levels=20)
# cbar = plt.colorbar(format='%.1e')
# cbar.set_label(r'Dissipation ($D^2 \mu^2 / \eta$)', rotation=270, labelpad=10)
# plt.xlabel(r'Horizontal   $x$ ($D$)', fontsize = 14) 
# plt.ylabel(r'  $y$ ($D$)', fontsize = 10) 
plt.title(r'(b) $T = 320\tau_M$, $\lambda = 2.5D$',loc ='left', fontsize = fs)
plt.yticks([])
plt.xticks([])
cbar = plt.colorbar(format='%.1e')
 
# cbar.set_label(r'Dissipation ($D^2 \mu^2 / \eta$)', rotation=270, labelpad=10)

  #####


  #####
w = wav[4]
TOPO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.topo_s_'
STRESS_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.stress.'
VELO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/velo'

TIME_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.time'
time = pd.read_table(TIME_FILEPATH, sep="\s+", skiprows=1, header = None)
time.columns = ['Timestep','Time', 'Total Dissipation','Total Elastic',
                    'Total Elastic+Dissipation','Total Work', 'Dissipation', 
                    'Elastic Energy','Work']

diss = time['Dissipation']
max_diss = max(diss)
ts = time.loc[time['Dissipation'] == max_diss, 'Timestep'].iloc[0]

if len(str(ts)) == 5:  t = str(ts)
elif len(str(ts)) == 4:  t = '0'+str(ts)
elif len(str(ts)) == 3:  t = '00'+str(ts)
elif len(str(ts)) == 2:  t = '000'+str(ts)
elif len(str(ts)) == 1:  t = '0000'+str(ts)



# read in relavent files for each time step i
topo = pd.read_table(TOPO_FILEPATH+str(t)+'.dat', sep="\s+", skiprows=1, header = None)
topo.columns = ['x','z','zprime','na','na']

stress = pd.read_table(STRESS_FILEPATH+str(ts), sep="\s+", skiprows=1, header = None)
stress.columns = ['SIGMA','Dissipation','ETA']
stress =stress[49:]
l = len(stress)

# visc = pd.read_table(VISC_FILEPATH, sep="\s+", skiprows=1, header = None)
# visc.columns = ['na','na','visc']

velo = pd.read_table(VELO_FILEPATH, sep="\s+", skiprows=1, header = None)
velo.columns = ['Node','X','Z','deltaX','deltaZ','na','na']
velo = velo[:l]  


x = velo['X']
y = velo['Z']
d = stress['Dissipation']/scale
# d[-1:] = 4.5e-10


    
plt.subplot2grid((3,3), (0,2)) 
plt.tricontourf(x,y-1,d, cmap=cmap,levels=20)
plt.yticks([])
plt.xticks([])
# plt.xlabel(r'Horizontal   $x$ ($D$)', fontsize = 14) 
# plt.ylabel(r'  $y$ ($D$)', fontsize = 10) 
plt.title(r'(c) $T = 320\tau_M$, $\lambda = 0.625D$',loc ='left', fontsize = fs)
cbar = plt.colorbar(format='%.1e')
cbar.set_label(r'dissipation ($\mu^2 / \eta$)', rotation=270, labelpad=15, fontsize=fs)
 
##############################
##############################

T = 'c05'
period = 20

wav = ['', 'A','B','C','D']

  #####
w = wav[2]
TOPO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.topo_s_'
STRESS_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.stress.'
VELO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/velo'

TIME_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.time'
time = pd.read_table(TIME_FILEPATH, sep="\s+", skiprows=1, header = None)
time.columns = ['Timestep','Time', 'Total Dissipation','Total Elastic',
                    'Total Elastic+Dissipation','Total Work', 'Dissipation', 
                    'Elastic Energy','Work']

diss = time['Dissipation']
max_diss = max(diss)
ts = time.loc[time['Dissipation'] == max_diss, 'Timestep'].iloc[0] 

if len(str(ts)) == 5:  t = str(ts)
elif len(str(ts)) == 4:  t = '0'+str(ts)
elif len(str(ts)) == 3:  t = '00'+str(ts)
elif len(str(ts)) == 2:  t = '000'+str(ts)
elif len(str(ts)) == 1:  t = '0000'+str(ts)



# read in relavent files for each time step i
topo = pd.read_table(TOPO_FILEPATH+str(t)+'.dat', sep="\s+", skiprows=1, header = None)
topo.columns = ['x','z','zprime','na','na']

stress = pd.read_table(STRESS_FILEPATH+str(ts), sep="\s+", skiprows=1, header = None)
stress.columns = ['SIGMA','Dissipation','ETA']
stress =stress[49:]
l = len(stress)

# visc = pd.read_table(VISC_FILEPATH, sep="\s+", skiprows=1, header = None)
# visc.columns = ['na','na','visc']

velo = pd.read_table(VELO_FILEPATH, sep="\s+", skiprows=1, header = None)
velo.columns = ['Node','X','Z','deltaX','deltaZ','na','na']
velo = velo[:l]  


x = velo['X']
y = velo['Z']
d = stress['Dissipation']/scale
# d[-1:] = 4.5e-10


    
plt.subplot2grid((3,3), (1,0))   
plt.tricontourf(x,y-1,d,cmap=cmap, levels=20)
plt.xticks([])
# cbar = plt.colorbar(format='%.1e')
# cbar.set_label(r'Dissipation ($D^2 \mu^2 / \eta$)', rotation=270, labelpad=10)
# plt.xlabel(r'Horizontal   $x$ ($D$)', fontsize = 14) 
plt.ylabel(r'  $y$ ($D$)', fontsize = fs) 
plt.title(r'(d) $T = 20\tau_M$, $\lambda = 10D$',loc ='left', fontsize = fs)
cbar = plt.colorbar(format='%.1e')
 
  #####



#####
w = wav[1]
TOPO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.topo_s_'
STRESS_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.stress.'
VELO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/velo'

TIME_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.time'
time = pd.read_table(TIME_FILEPATH, sep="\s+", skiprows=1, header = None)
time.columns = ['Timestep','Time', 'Total Dissipation','Total Elastic',
                    'Total Elastic+Dissipation','Total Work', 'Dissipation', 
                    'Elastic Energy','Work']

diss = time['Dissipation']
max_diss = max(diss)
ts = time.loc[time['Dissipation'] == max_diss, 'Timestep'].iloc[0] 

if len(str(ts)) == 5:  t = str(ts)
elif len(str(ts)) == 4:  t = '0'+str(ts)
elif len(str(ts)) == 3:  t = '00'+str(ts)
elif len(str(ts)) == 2:  t = '000'+str(ts)
elif len(str(ts)) == 1:  t = '0000'+str(ts)



# read in relavent files for each time step i
topo = pd.read_table(TOPO_FILEPATH+str(t)+'.dat', sep="\s+", skiprows=1, header = None)
topo.columns = ['x','z','zprime','na','na']

stress = pd.read_table(STRESS_FILEPATH+str(ts), sep="\s+", skiprows=1, header = None)
stress.columns = ['SIGMA','Dissipation','ETA']
stress =stress[49:]
l = len(stress)

# visc = pd.read_table(VISC_FILEPATH, sep="\s+", skiprows=1, header = None)
# visc.columns = ['na','na','visc']

velo = pd.read_table(VELO_FILEPATH, sep="\s+", skiprows=1, header = None)
velo.columns = ['Node','X','Z','deltaX','deltaZ','na','na']
velo = velo[:l]  


x = velo['X']
y = velo['Z']
d = stress['Dissipation']/scale
# d[-1:] = 4.5e-10


    
plt.subplot2grid((3,3), (1,1)) 
plt.tricontourf(x,y-1,d, cmap=cmap,levels=20)
plt.xticks([])
# cbar = plt.colorbar(format='%.1e')
# cbar.set_label(r'Dissipation ($D^2 \mu^2 / \eta$)', rotation=270, labelpad=10)
# plt.xlabel(r'Horizontal   $x$ ($D$)', fontsize = 14) 
# plt.ylabel(r'  $y$ ($D$)', fontsize = 10) 
plt.title(r'(e) $T = 20\tau_M$, $\lambda = 2.5D$',loc ='left', fontsize = fs)
cbar = plt.colorbar(format='%.1e')
plt.yticks([])
 
# cbar.set_label(r'Dissipation ($D^2 \mu^2 / \eta$)', rotation=270, labelpad=10)

  #####

  #####
w = wav[4]
TOPO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.topo_s_'
STRESS_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.stress.'
VELO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/velo'

TIME_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.time'
time = pd.read_table(TIME_FILEPATH, sep="\s+", skiprows=1, header = None)
time.columns = ['Timestep','Time', 'Total Dissipation','Total Elastic',
                    'Total Elastic+Dissipation','Total Work', 'Dissipation', 
                    'Elastic Energy','Work']

diss = time['Dissipation']
max_diss = max(diss)
ts = time.loc[time['Dissipation'] == max_diss, 'Timestep'].iloc[0] -90

if len(str(ts)) == 5:  t = str(ts)
elif len(str(ts)) == 4:  t = '0'+str(ts)
elif len(str(ts)) == 3:  t = '00'+str(ts)
elif len(str(ts)) == 2:  t = '000'+str(ts)
elif len(str(ts)) == 1:  t = '0000'+str(ts)



# read in relavent files for each time step i
topo = pd.read_table(TOPO_FILEPATH+str(t)+'.dat', sep="\s+", skiprows=1, header = None)
topo.columns = ['x','z','zprime','na','na']

stress = pd.read_table(STRESS_FILEPATH+str(ts), sep="\s+", skiprows=1, header = None)
stress.columns = ['SIGMA','Dissipation','ETA']
stress =stress[49:]
l = len(stress)

# visc = pd.read_table(VISC_FILEPATH, sep="\s+", skiprows=1, header = None)
# visc.columns = ['na','na','visc']

velo = pd.read_table(VELO_FILEPATH, sep="\s+", skiprows=1, header = None)
velo.columns = ['Node','X','Z','deltaX','deltaZ','na','na']
velo = velo[:l]  


x = velo['X']
y = velo['Z']
d = stress['Dissipation']/scale
# d[-1:] = 4.5e-10


    
plt.subplot2grid((3,3), (1,2)) 
plt.tricontourf(x,y-1,d,cmap=cmap, levels=20)
 
plt.yticks([])
plt.xticks([])
# plt.xlabel(r'Horizontal   $x$ ($D$)', fontsize = 14) 
# plt.ylabel(r'  $y$ ($D$)', fontsize = 10) 
plt.title(r'(f) $T = 20\tau_M$, $\lambda = 0.625D$',loc ='left', fontsize = fs)
cbar = plt.colorbar(format='%.1e')
cbar.set_label(r'dissipation ($\mu^2 / \eta$)', rotation=270, labelpad=15, fontsize=fs)

#############################
############################

T = 'c15'
period = 320
ts = 156
t = '00156'
wav = ['', 'A','B','C','D']

  #####
w = wav[2]
TOPO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.topo_s_'
STRESS_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.stress.'
VELO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/velo'

TIME_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.time'
time = pd.read_table(TIME_FILEPATH, sep="\s+", skiprows=1, header = None)
time.columns = ['Timestep','Time', 'Total Dissipation','Total Elastic',
                    'Total Elastic+Dissipation','Total Work', 'Dissipation', 
                    'Elastic Energy','Work']

diss = time['Dissipation']
max_diss = max(diss)
ts = time.loc[time['Dissipation'] == max_diss, 'Timestep'].iloc[0] 

if len(str(ts)) == 5:  t = str(ts)
elif len(str(ts)) == 4:  t = '0'+str(ts)
elif len(str(ts)) == 3:  t = '00'+str(ts)
elif len(str(ts)) == 2:  t = '000'+str(ts)
elif len(str(ts)) == 1:  t = '0000'+str(ts)



# read in relavent files for each time step i
topo = pd.read_table(TOPO_FILEPATH+str(t)+'.dat', sep="\s+", skiprows=1, header = None)
topo.columns = ['x','z','zprime','na','na']

stress = pd.read_table(STRESS_FILEPATH+str(ts), sep="\s+", skiprows=1, header = None)
stress.columns = ['SIGMA','Dissipation','ETA']
stress =stress[49:]
l = len(stress)

# visc = pd.read_table(VISC_FILEPATH, sep="\s+", skiprows=1, header = None)
# visc.columns = ['na','na','visc']

velo = pd.read_table(VELO_FILEPATH, sep="\s+", skiprows=1, header = None)
velo.columns = ['Node','X','Z','deltaX','deltaZ','na','na']
velo = velo[:l]  


x = velo['X']
y = velo['Z']
d = stress['Dissipation']/scale
# d[-1:] = 4.5e-10


    
plt.subplot2grid((3,3), (2,0)) 
plt.tricontourf(x,y-1,d,cmap=cmap, levels=20)
# cbar = plt.colorbar(format='%.1e')
# cbar.set_label(r'Dissipation ($D^2 \mu^2 / \eta$)', rotation=270, labelpad=10)
plt.xlabel(r'  $x$ ($D$)', fontsize = fs) 
plt.ylabel(r'  $y$ ($D$)', fontsize = fs) 
plt.title(r'(g) $T = 0.02 \tau_M$, $\lambda = 10D$',loc ='left', fontsize = fs)
cbar = plt.colorbar(format='%.1e')
 
  #####


  #####
w = wav[1]
TOPO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.topo_s_'
STRESS_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.stress.'
VELO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/velo'

TIME_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.time'
time = pd.read_table(TIME_FILEPATH, sep="\s+", skiprows=1, header = None)
time.columns = ['Timestep','Time', 'Total Dissipation','Total Elastic',
                    'Total Elastic+Dissipation','Total Work', 'Dissipation', 
                    'Elastic Energy','Work']

diss = time['Dissipation']
max_diss = max(diss)
ts = time.loc[time['Dissipation'] == max_diss, 'Timestep'].iloc[0]+6

if len(str(ts)) == 5:  t = str(ts)
elif len(str(ts)) == 4:  t = '0'+str(ts)
elif len(str(ts)) == 3:  t = '00'+str(ts)
elif len(str(ts)) == 2:  t = '000'+str(ts)
elif len(str(ts)) == 1:  t = '0000'+str(ts)


# read in relavent files for each time step i
topo = pd.read_table(TOPO_FILEPATH+str(t)+'.dat', sep="\s+", skiprows=1, header = None)
topo.columns = ['x','z','zprime','na','na']

stress = pd.read_table(STRESS_FILEPATH+str(ts), sep="\s+", skiprows=1, header = None)
stress.columns = ['SIGMA','Dissipation','ETA']
stress =stress[49:]
l = len(stress)

# visc = pd.read_table(VISC_FILEPATH, sep="\s+", skiprows=1, header = None)
# visc.columns = ['na','na','visc']

velo = pd.read_table(VELO_FILEPATH, sep="\s+", skiprows=1, header = None)
velo.columns = ['Node','X','Z','deltaX','deltaZ','na','na']
velo = velo[:l]  


x = velo['X']
y = velo['Z']
d = stress['Dissipation']/scale
# d[-1:] = 4.5e-10


    
plt.subplot2grid((3,3), (2,1))  
plt.tricontourf(x,y-1,d,cmap=cmap, levels=20)
plt.yticks([])
# cbar = plt.colorbar(format='%.1e')
# cbar.set_label(r'Dissipation ($D^2 \mu^2 / \eta$)', rotation=270, labelpad=10)
plt.xlabel(r'  $x$ ($D$)', fontsize = fs) 
# plt.ylabel(r'  $y$ ($D$)', fontsize = 10) 
plt.title(r'(h) $T = 0.02 \tau_M$, $\lambda = 2.5D$',loc ='left', fontsize = fs)
cbar = plt.colorbar(format='%.1e')
# cbar.set_label(r'Dissipation ($D^2 \mu^2 / \eta$)', rotation=270, labelpad=10)
 




  #####
w = wav[4]
TOPO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.topo_s_'
STRESS_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.stress.'
VELO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/velo'

TIME_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.time'
time = pd.read_table(TIME_FILEPATH, sep="\s+", skiprows=1, header = None)
time.columns = ['Timestep','Time', 'Total Dissipation','Total Elastic',
                    'Total Elastic+Dissipation','Total Work', 'Dissipation', 
                    'Elastic Energy','Work']

diss = time['Dissipation']
max_diss = max(diss)
ts = time.loc[time['Dissipation'] == max_diss, 'Timestep'].iloc[0] +6

if len(str(ts)) == 5:  t = str(ts)
elif len(str(ts)) == 4:  t = '0'+str(ts)
elif len(str(ts)) == 3:  t = '00'+str(ts)
elif len(str(ts)) == 2:  t = '000'+str(ts)
elif len(str(ts)) == 1:  t = '0000'+str(ts)



# read in relavent files for each time step i
topo = pd.read_table(TOPO_FILEPATH+str(t)+'.dat', sep="\s+", skiprows=1, header = None)
topo.columns = ['x','z','zprime','na','na']

stress = pd.read_table(STRESS_FILEPATH+str(ts), sep="\s+", skiprows=1, header = None)
stress.columns = ['SIGMA','Dissipation','ETA']
stress =stress[49:]
l = len(stress)

# visc = pd.read_table(VISC_FILEPATH, sep="\s+", skiprows=1, header = None)
# visc.columns = ['na','na','visc']

velo = pd.read_table(VELO_FILEPATH, sep="\s+", skiprows=1, header = None)
velo.columns = ['Node','X','Z','deltaX','deltaZ','na','na']
velo = velo[:l]  


x = velo['X']
y = velo['Z']
d = stress['Dissipation']/scale
# d[-1:] = 4.5e-10


    
plt.subplot2grid((3,3), (2,2)) 
plt.tricontourf(x,y-1,d,cmap=cmap, levels=20)
plt.yticks([])
plt.xticks([0.0,0.15,0.3])
plt.xlabel(r'  $x$ ($D$)', fontsize = fs) 
# plt.ylabel(r'  $y$ ($D$)', fontsize = 10) 
plt.title(r'(i) $T = 0.02 \tau_M$, $\lambda = 0.625D$',loc ='left', fontsize = fs)
cbar = plt.colorbar(format='%.1e')
cbar.set_label(r'dissipation ($\mu^2 / \eta$)', rotation=270, labelpad=15, fontsize=fs)
 


plt.tight_layout()
plt.savefig('/Users/emmadevin/School/ViscousDissipation/Paper_Figures/Figure_8.eps')

