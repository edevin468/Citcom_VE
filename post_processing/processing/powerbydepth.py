#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Emma Devin
This program takes dissipation data for every node and time averages the power for every row.
"""

import numpy as np
import pandas as pd
float_formatter = '{:.1f}'.format

# wavelength labels and values
wav = ['','A','B','C','D']

wavelengths = {'':2.5, 'A':1.25,'B':5.0,'C':0.625,'D':0.3125} 

# List of viscosity structures
cases = ['100']#,'500']           


# number of nodes in horizontal direction by wavelength
sizes = {'':97, 'A':65,'B':193,'C':33,'D':33}


# list of periods
periods = [0,320,160,80,40,20,10,5,2.5,1.25,0.625,0.3125,0.15625,0.078125,0.0390625,
            0.0195,0.00976,0.00488,0.00244,0.00122,0.000610,0.000305,0.000152,0.0000763]


# list of number of timesteps for each period
lengths = [0,612,612,3030,3030,1515,612,612,612,312,312,312,312,312,312,312,312,312,
            312,312,312,312,312,312]

# loop through all cases and wavelengths and periods             
for c in cases: 
    for w in wav:
        maxdf = pd.DataFrame()
        df1 = pd.DataFrame() 
        
        # read relevant files into Python
        VELO_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+c+w+'_P/velo'
        
        velo = pd.read_table(VELO_FILEPATH, sep="\s+", skiprows=1, header = None)
        velo.columns = ['Node','X','Z','deltaX','deltaZ','na','na']
        velo = velo[:49]
        z = velo['Z']
        
    
    
        for p in range(23):
        
            # define string for period label
            p+=1
            if len(str(p))==1: T = 'c0'+str(p)
            elif len(str(p))==2: T = 'c'+str(p)
            
            # define step size based on period and number of steps
            steps = lengths[p]
            if steps == 312: stepsize = 12
            elif steps == 612: stepsize = 12
            elif steps ==3030: stepsize = 60
            elif steps == 1515: stepsize = 30
            
            period = periods[p]
     
            df = pd.DataFrame()
            df2 = pd.DataFrame()
            for i in range(0,steps,stepsize):
                
                t = i*period/steps 
         
                STRESS_FILEPATH = '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+c+w+'_P/'+T+'.stress.'
                    
                # define index j for reading in topo files for each timestep
                if len(str(i))==1: j = '0000'+str(i)
                elif len(str(i))==3: j = '00'+str(i) 
                elif len(str(i))==4: j = '0'+str(i)
                
                stress = pd.read_table(STRESS_FILEPATH+str(i), sep="\s+", skiprows=1, header = None)
                stress.columns = ['SIGMA','Dissipation','ETA']
                stress =stress[49:]
                l = len(stress)
                
                d = stress['Dissipation']
            
                s = sizes[w]
                
                d = np.array(d) 
                d = d.reshape(s, 49) 
                d = np.transpose(d)
                d = pd.DataFrame(d)
                
                # did the math! for real this time -- reduces down to just the average, or equivalently this
                # look at photo of math before messing with this again
                total_diss = d.sum(axis=1)/sizes[w]
                df[str(t)] = total_diss
            
                df2[i] = total_diss
            
            df_t = df.T
            df_t['time']=df_t.index
            df_t = df_t[1:]
            
            LAYER_POWER = []
            
            # compute time averaged power for each row of nodes
            for i in range(49):
     
                layer = df_t[i]
                time = df_t.index.tolist()
                time = time[1:]
                t = [float(i) for i in time] 
                l = len(t)
          
                j = 1
                d_layer = 0
     
                while j<l:
                    layer_curr = layer[j]
                    layer_prev = layer[j-1]
                    time_curr = t[j]
                    time_prev = t[j-1] 
                    
                    d_layer += 0.5*(layer_curr + layer_prev)*(time_curr - time_prev)
                
                    j += 1
             
                layer_pow_avg = d_layer/(period)
                LAYER_POWER.append(layer_pow_avg)
    
            df1[str(period)] = LAYER_POWER
            
                                    
        df1.to_csv('/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power_by_depth'
                    '/Avg_'+c+'/Case'+c+w+'.csv')
        