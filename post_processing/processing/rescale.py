#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Written by Emma Devin
This program takes data from files containing spatially averaged dissipation, work, and elastic 
energy data for each case and each timestep and rescales them to real units for several real processes
'''

import numpy as np
import pandas as pd

viscosities = [0,1e22,1e21,1e20,1e19,1e18]        # list of reference viscosities 
# viscosities = [1e18]

# VS1 = '100', VS2 = '200', and VS3 = '500'
cases = ['900']                                   # designates viscosity structure
wav = ['', 'A','B','C','D']                     # list of wavelength labels 

# wavelength labels and their values
wavelengths = {'':2.5, 'A':1.25,'B':5.0,'C':0.625,'D':0.3125}

for c in cases:

    v = 1
    while v<=5:
        eta = viscosities[v]                        # reference viscosity in Pa s
        mu = 7e10                                   # reference shear modulus in Pa                            
        D = 2.9e6                                   # box depth in m 
          
        power_scale = mu**2*D/eta                # scaling for power in W/m^2
        loading_scale = 290                         # scaling for 1m load
        ice_load = 1000                             # scaling for 1000m rock load to simulate 3000m ice
        tide_load = 1000
        amazon_load = 0.3                           # scaling for 0.3m rock load to simulate 30cm water
        ocean_load = 1                              # scaling for OCEAN ***********
        maxwell_time = eta/mu                       # maxwell time in s
    
        for w in wav:
        
            # make string sequence to read files for each loading period
            k = 1
            while k < 24:
                if len(str(k))== 1:
                    i = 'c0'+str(k)
                else:
                    i = 'c'+str(k)
                k+=1
                
                # read files into Python
                FILE =  '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+c+w+'/'+i+'.time'          
                
                df = pd.read_table(FILE, sep="\s+", skiprows=1, header = None)
                df.columns = ['Timestep','Time', 'Total Dissipation','Total Elastic',
                    'Total Elastic+Dissipation','Total Work', 'Dissipation', 
                    'Elastic Energy','Work']
                
                # output in *.time files is dissipation/elastic energy integrated over area 
                # and work integrated over box width so we must divide out those quantities
                # to get power 
                
                area = 1*wavelengths[w]
                
                t = np.array(df['Time'])
                dissipation = df['Dissipation']/area
                elastic = df['Elastic Energy']/area
                work = df['Work']/wavelengths[w]
                
                # loading scales are squared because energy goes like the 
                # square of the displacement
                # all of these outputs will be in W/m^2
                
                # used the next three lines for 1m rock load
                real_diss = (dissipation*power_scale)/loading_scale**2         
                real_work = (work*power_scale)/loading_scale**2                
                real_elastic = (elastic*power_scale)/loading_scale**2
                
                # real_diss = (dissipation*power_scale*tide_load**2)/loading_scale**2         
                # real_work = (work*power_scale*tide_load**2)/loading_scale**2                
                # real_elastic = (elastic*power_scale*tide_load**2)/loading_scale**2
        
                # use the next three lines for ice load
                # real_diss = (dissipation*power_scale*ice_load**2)/loading_scale**2         
                # real_work = (work*power_scale*ice_load**2)/loading_scale**2                
                # real_elastic = (elastic*power_scale*ice_load**2)/loading_scale**2    
                
                # use the next three lines for river basin load
                # real_diss = (dissipation*power_scale*amazon_load**2)/loading_scale**2        
                # real_work = (work*power_scale*amazon_load**2)/loading_scale**2              
                # real_elastic = (elastic*power_scale*amazon_load**2)/loading_scale**2    
                
                # use the next three lines for ocean load
                # real_diss = (dissipation*power_scale*ocean_load**2)/loading_scale**2        
                # real_work = (work*power_scale*ocean_load**2)/loading_scale**2              
                # real_elastic = (elastic*power_scale*ocean_load**2)/loading_scale**2   
                
                real_time = t*maxwell_time                    
                
                # save data to *.csv files
                df1 = pd.DataFrame()
                df1['Time'] = real_time
                df1['Dissipation'] = real_diss
                df1['Work'] = real_work
                df1['Elastic'] = real_elastic
                
                # use this filepath for 1m rock load
                df1.to_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/'
                    'Viscosity_'+str(v)+'/Total_Dissipation/Case_'+c+'/'
                    'CASE'+c+w+i+'.csv')
                
                # df1.to_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/'
                #     'Viscosity_'+str(v)+'/Total_Dissipation_tides/Case_'+c+'/'
                #     'CASE'+c+w+i+'.csv')
                
                # use this filepath for ice load
                # df1.to_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/'
                #     'Viscosity_'+str(v)+'/Total_Dissipation_ice/Case_'+c+'/'
                #     'CASE'+c+w+i+'.csv')
                
                # # use this filepath for river basin load    
                # df1.to_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/'
                #     'Viscosity_'+str(v)+'/Total_Dissipation_amazon/Case_'+c+'/'
                #     'CASE'+c+w+i+'.csv')
                
                # # use this filepath for ocean load    
                # df1.to_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/'
                #     'Viscosity_'+str(v)+'/Total_Dissipation_ocean/Case_'+c+'/'
                #     'CASE'+c+w+i+'.csv')
        v+=1