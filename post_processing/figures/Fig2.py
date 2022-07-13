#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 16:18:01 2021

@author: emmadevin
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

fig = plt.figure(constrained_layout=True, figsize=(10,8))
gs = fig.add_gridspec(nrows=2,ncols=6)
ax1=fig.add_subplot(gs[0,0:2])
ax2=fig.add_subplot(gs[0,2:4])
ax3=fig.add_subplot(gs[0,4:6])
ax4=fig.add_subplot(gs[1,0:3])
ax5=fig.add_subplot(gs[1,3:6])

fig.patch.set_facecolor('white')
plt.style.use('classic')


case = '100'
wav = ['', 'A','B','C','D']
wavelengths = {'':2.5, 'A':1.25,'B':5.0,'C':0.625,'D':0.3125}

periods = [0,320,160,80,40,20,10,5,2.5,1.25,0.625,0.3125,0.15625,0.078125,0.0390625,0.0195,0.00976,0.00488,0.00244,0.00122,0.000610,0.000305,0.000152,0.0000763]

# 
w = 'C'
sigma_0 = 3300*10*1
T = 'c04'
k = 4
s = 600
period = periods[k]
    
if k < 3:
    m = 612
    s = 12
elif k >= 3 and k < 5:
    m = 3020
    s = 30
elif k >= 5 and k < 6:
    m = 1515
    s = 15
elif k >= 6 and k < 9: 
    m = 612
    s = 12
else:
    m = 312
    s = 12


#=====================#
# panel (a)
#=====================#

i = 0
t = 0

x1 = np.linspace(0,wavelengths[w],300)
y1 = sigma_0*np.cos(2*np.pi*x1/1.25)*np.sin(2*np.pi*t/40)


TOPO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'_P/'+T+'.topo_s_'

# define index j for reading in topo files
if len(str(i))==1: j = '0000'+str(i)
elif len(str(i))==2: j = '000'+str(i) 
elif len(str(i))==3: j = '00'+str(i) 
elif len(str(i))==4: j = '0'+str(i)

k = i + 30
if len(str(k))==1: m = '0000'+str(k)
elif len(str(k))==2: m = '000'+str(k) 
elif len(str(k))==3: m = '00'+str(k) 
elif len(str(k))==4: m = '0'+str(k)

# read in relavent files for each time step i
topo = pd.read_table(TOPO_FILEPATH+str(j)+'.dat', sep="\s+", skiprows=0, header = None)
topo.columns = ['x','z','zprime','na','na']

topo1 = pd.read_table(TOPO_FILEPATH+str(m)+'.dat', sep="\s+", skiprows=0, header = None)
topo1.columns = ['x','z','zprime','na','na']


X = topo['x']
z = topo['z']*1e4
z1 = topo1['z']*1e4
disp = z1-z

ax1.set_xlabel(r'x ($D$)')
ax1.set_ylabel('topography (m)')
ax1.plot(X,z, c = 'k',ls='--',lw=1, label='topography')
ax1.plot(X,disp,c = 'k', ls=':',lw=2,label='incremental displacement')

ax1.set_xlim(0,0.625)
ax1.set_ylim(-1,1)
ax1.set_xticks([0.0,0.3,0.6])

ax1b = ax1.twinx()
# ax1b.set_ylabel('loading force (Pa)')
ax1b.plot(x1,y1/1000, c = 'k', ls='-',lw=1, label = 'loading force ')
ax1b.plot(-x1,y1, c = 'k', ls='--',lw=1, label= 'topography') # plot off screen to get info for legend
ax1b.plot(-X,disp,c = 'k', ls=':',lw=2, label='incremental \ndisplacement') # plot off screen to get info for legend
ax1b.set_ylim(-40,40)
ax1b.set_yticklabels([])

plt.title(r'(a) t = 0$\tau_M$',loc ='left', fontsize =12, color = 'k',pad=15 )
ax1b.legend(fontsize=10)

#=====================#
# panel (b)
#=====================#

# i = 60
i = 360

t = 5

x1 = np.linspace(0,wavelengths[w],300)
y1 = sigma_0*np.cos(2*np.pi*x1/1.25)*np.sin(2*np.pi*t/40)


TOPO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'_P/'+T+'.topo_s_'

# define index j for reading in topo files
if len(str(i))==1: j = '0000'+str(i)
elif len(str(i))==2: j = '000'+str(i) 
elif len(str(i))==3: j = '00'+str(i) 
elif len(str(i))==4: j = '0'+str(i)


k = i + 30
if len(str(k))==1: m = '0000'+str(k)
elif len(str(k))==2: m = '000'+str(k) 
elif len(str(k))==3: m = '00'+str(k) 
elif len(str(k))==4: m = '0'+str(k)


# read in relavent files for each time step i
topo = pd.read_table(TOPO_FILEPATH+str(j)+'.dat', sep="\s+", skiprows=0, header = None)
topo.columns = ['x','z','zprime','na','na']

topo1 = pd.read_table(TOPO_FILEPATH+str(m)+'.dat', sep="\s+", skiprows=0, header = None)
topo1.columns = ['x','z','zprime','na','na']


X = topo['x']
z = topo['z']*1e4
z1 = topo1['z']*1e4
disp = z1-z


ax2.set_xlabel(r'x ($D$)')
# ax2.set_ylabel('topography (m)')
ax2.plot(X,z, c = 'k',ls='--',lw=1,label='topography')
ax2.plot(X,disp,c = 'k', ls=':',lw=2, label='incremental displacement')
ax2.set_xlim(0,0.625)
ax2.set_ylim(-1,1)
ax2.set_xticks([0.0,0.3,0.6])
ax2.set_yticklabels([])


ax2b = ax2.twinx()
# ax2b.set_ylabel('loading force (Pa)')
ax2b.plot(x1,y1/1000, c = 'k',ls='-',lw=1, label = 'loading force ')

ax2b.set_ylim(-40,40)
ax2b.set_yticklabels([])

ax2.set_title(r'(b) t = 5$\tau_M$',loc ='left', fontsize =12, color = 'k',pad=15 )


# for text in leg.get_texts():
#     plt.setp(text, color = 'k')


#=====================#
# panel (c)
#=====================#

# i = 168
i = 750

t = 10

x1 = np.linspace(0,wavelengths[w],300)
y1 = sigma_0*np.cos(2*np.pi*x1/1.25)*np.sin(2*np.pi*t/40)


TOPO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'_P/'+T+'.topo_s_'

# define index j for reading in topo files
if len(str(i))==1: j = '0000'+str(i)
elif len(str(i))==2: j = '000'+str(i) 
elif len(str(i))==3: j = '00'+str(i) 
elif len(str(i))==4: j = '0'+str(i)


k = i + 30
if len(str(k))==1: m = '0000'+str(k)
elif len(str(k))==2: m = '000'+str(k) 
elif len(str(k))==3: m = '00'+str(k) 
elif len(str(k))==4: m = '0'+str(k)


# read in relavent files for each time step i
topo = pd.read_table(TOPO_FILEPATH+str(j)+'.dat', sep="\s+", skiprows=0, header = None)
topo.columns = ['x','z','zprime','na','na']

topo1 = pd.read_table(TOPO_FILEPATH+str(m)+'.dat', sep="\s+", skiprows=0, header = None)
topo1.columns = ['x','z','zprime','na','na']


X = topo['x']
z = topo['z']*1e4
z1 = topo1['z']*1e4
disp = z1-z


ax3.set_xlabel(r'x ($D$)')
# ax3.set_ylabel('topography (m)')
ax3.plot(X,z, c = 'k',ls='--',lw=1, label='topography')
ax3.plot(X,disp,c = 'k',ls=':',lw=2, label='incremental displacement')
ax3.set_xlim(0,0.625)
ax3.set_ylim(-1,1)
ax3.set_xticks([0.0,0.3,0.6])
ax3.set_yticklabels([])

ax3b = ax3.twinx()
ax3b.set_ylabel('loading force (KPa)',rotation= 270,labelpad=15)
ax3b.plot(x1,y1/1000, c = 'k',ls='-',lw=1, label = 'loading force ')


ax3b.set_ylim(-40,40)
# plt.yticks(rotation=270)

ax3.set_title(r'(c) t = 10$\tau_M$',loc ='left', fontsize =12, color = 'k',pad=15 )

# ax3b.legend(loc='right', bbox_to_anchor=(2.5, 0.8), ncol = 1, fontsize = 10, facecolor= 'w')

#=====================#
# panel (d)
#=====================#

case = '100'
wav = ['', 'A','B','C','D']
wavelengths = {'':2.5, 'A':1.25,'B':5.0,'C':0.625,'D':0.3125}


periods = [0,320,160,80,40,20,10,5,2.5,1.25,0.625,0.3125,0.15625,0.078125,0.0390625,0.0195,0.00976,0.00488,0.00244,0.00122,0.000610,0.000305,0.000152,0.0000763]

max_topo = []
max_load = []
time = []
disp = []

w = 'C'
area = 1*wavelengths[w]
h =1e-4
T = 'c04'
k = 4
period = periods[k]
    
if k < 3:
    m = 612
    s = 12
elif k >= 3 and k < 5:
    m = 3020
    s = 30
elif k >= 5 and k < 6:
    m = 1515
    s = 15
elif k >= 6 and k < 9: 
    m = 612
    s = 12
else:
    m = 312
    s = 12
    

for i in range(0,m,s):
    
    t = (i)*period/(m-s)
    
    x1 = np.linspace(0,wavelengths[w],600)

    
    TOPO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'_P/'+T+'.topo_s_'
    # TOPO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+T+'.topo_s_'
    
    # define index j for reading in topo files
    if len(str(i))==1: j = '0000'+str(i)
    elif len(str(i))==2: j = '000'+str(i) 
    elif len(str(i))==3: j = '00'+str(i) 
    elif len(str(i))==4: j = '0'+str(i)
    
    # read in relavent files for each time step i
    topo = pd.read_table(TOPO_FILEPATH+str(j)+'.dat', sep="\s+", skiprows=0, header = None)
    topo.columns = ['x','z','zprime','na','na']
    
    
    X = topo['x']
    z = topo['z']*1e4
    
    max_z = z[0]

    max_topo.append(max_z)

    time.append(t)

    
max_topo = np.array(max_topo)

for j in range(len(max_topo)-1):
    inc_displacement = max_topo[j+1]-max_topo[j]
    disp.append(inc_displacement)

disp.insert(0,0)

ax4.plot(time, max_topo, c = 'k',lw=1,ls='--',label='topography')
ax4.plot(time, disp, c='k',lw=2,ls=':',label = 'incremental \ndisplacement')
ax4.set_xlabel(r'time ($\tau_M$)', fontsize =11)
ax4.set_ylabel(r'y (m)', fontsize = 11)
ax4.set_xlim(0,40)
ax4.set_ylim(-0.8,0.8)
ax4.set_title(r'(d)' ,loc ='left', fontsize =12, color = 'k',pad=15 )
# ax4.set_xticks([0,,40])
ax4.legend(fontsize=10,loc='upper left')



#=====================#
# panel (e)
#=====================#


case = '100'

i = 'c04'

wav = ['', 'A','B','C','D']
box_widths = {'':2.5, 'A':1.25,'B':5.0,'C':0.625,'D':0.3125}
w = 'C'

FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Total_Dissipation/Case_'+case+'/CASE'+case+w+i+'.csv'

df = pd.read_csv(FILE)
l = len(df)

box_width = box_widths[w]       
t = np.array(df['Time'])
# t = [float(float_formatter(x)) for x in t]
dissipation = df['Dissipation']*1
elastic = df['Elastic']*1
work = df['Work']*1
balance = dissipation + elastic + work

scale = 290**2 # scale for dimensionless 1m rock load

diss = dissipation/scale
el = elastic/scale
w = work/scale
bal = balance/scale

T = max(t) 


ax5.plot(t, bal, c ='k', ls = '--',lw = 2,label = 'sum')
ax5.plot(t, diss, c ='k',ls = '-', lw = 2,label ='dissipation')
ax5.plot(t, el, c = 'k' ,ls = '-', lw = 1, label = 'elastic energy')
ax5.yaxis.tick_right()
ax5.yaxis.set_label_position("right")
ax5.plot(t, w, c='k'    ,ls = ':', lw = 2, label = 'work')
ax5.set_xlim(0,40)
ax5.set_ylim(-1.8*10**-14,1.8*10**-14)
ax5.set_xlabel(r'time ($\tau_M$)', fontsize = 11)
ax5.set_ylabel(r'energy flux ($D \mu^2 / \eta$)', fontsize = 11,rotation=270, labelpad=15)
ax5.set_title(r'(e)' ,loc ='left', fontsize =12, color = 'k',pad=15 )
ax5.legend(fontsize=10,loc='upper left', ncol=1)
# leg = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), facecolor = 'w', ncol = 2, fontsize = 10)

fig.savefig('/Users/emmadevin/School/ViscousDissipation/Paper_Figures/Figure_2.eps', bbox_inches = 'tight')
