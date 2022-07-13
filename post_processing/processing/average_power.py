#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Emma Devin
This program takes rescaled or dimensionless spatially averaged power data 
and calculates the time-averaged power for each period.
"""

import numpy as np
import pandas as pd

viscosities = [0,1e22,1e21,1e20,1e19,1e18]        # list of reference viscosities 

# VS1 = '100', VS2 = '200', and VS3 = '500'
cases = ['900']                      # designates viscosity structure
wav = ['', 'A','B','C','D']                     # list of wavelength labels 

# wavelength labels and their values
wavelengths = {'':2.5, 'A':1.25,'B':5.0,'C':0.625,'D':0.3125}

for c in cases: 
    v = 1
    while v<=5:
        for w in wav:
            DISS_POWER = []
            WORK_POWER = []
            ELASTIC_POWER = []
            PERIOD =[]
            
            k = 1
            while k < 24:
                if len(str(k))== 1:
                    i = 'c0'+str(k)
                else:
                    i = 'c'+str(k)
                k+=1
                
                # read data rescaled by rescale.py saved to .csv
                
                # use this filepath for dimensionless data
                # FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Total_Dissipation/Case_'+c+'/CASE'+c+w+i+'.csv'
                
                # use this filepath for 1m rock load
                FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_'+str(v)+'/Total_Dissipation/Case_'+c+'/CASE'+c+w+i+'.csv'
                
                # use this filepath for ice load
                # FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_'+str(v)+'/Total_Dissipation_ice/Case_'+c+'/CASE'+c+w+i+'.csv'
                
                # use this filepath for river basin load
                # FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_'+str(v)+'/Total_Dissipation_amazon/Case_'+c+'/CASE'+c+w+i+'.csv'
               
                # use this filepath for tide load
                # FILE = '/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_'+str(v)+'/Total_Dissipation_tides/Case_'+c+'/CASE'+c+w+i+'.csv'
    
                df = pd.read_csv(FILE)
                l = len(df)
                
                t = np.array(df['Time'])
                dissipation = df['Dissipation']
                elastic = df['Elastic']
                work = df['Work']
                period =  t[l-1]
            
                j = 1
                diss_pow = 0
                work_pow = 0
                elastic_pow = 0
                
                # calculate time integral using trapezoids
                while j<l:
                 
                    diss_curr = dissipation[j]
                    diss_prev = dissipation[j-1]
                
                    work_curr = work[j]
                    work_prev = work[j-1]
                          
                    elastic_curr = elastic[j]
                    elastic_prev = elastic[j-1]
                
                    time_curr = t[j]
                    time_prev = t[j-1] 
                 
                    diss_pow += 0.5*(diss_curr + diss_prev)*(time_curr - time_prev)
                    work_pow += 0.5*(work_curr + work_prev)*(time_curr - time_prev)
                    elastic_pow += 0.5*(elastic_curr + elastic_prev)*(time_curr - time_prev)
                    
                    j += 1
               
                # divided by period to get time average
                diss_pow_avg = diss_pow/period 
                work_pow_avg = work_pow/period
                elastic_pow_avg = elastic_pow/period
            
                DISS_POWER.append(diss_pow_avg)
                WORK_POWER.append(work_pow_avg)
                ELASTIC_POWER.append(elastic_pow_avg)
                    
                PERIOD.append(period)
          
            # save to relevant folder
            df1 = pd.DataFrame() 
            df1['Period'] = PERIOD
            df1['Dissipative_Power'] = DISS_POWER
            df1['Work_Power'] = WORK_POWER
            df1['ELastic_Power'] = ELASTIC_POWER
            
            # use this filepath for dimensionless data
            # df1.to_csv('/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Avg_Power/Case_'+c+'/Case'+c+w+'.csv')
           
            # use this filepath for 1m rock load
            df1.to_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_'+str(v)+'/Avg_Power/Case_'+c+'/Case'+c+w+'.csv') 
           
            # use this filepath for ice load
            # df1.to_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_'+str(v)+'/Avg_Power_ice/Case_'+c+'/Case'+c+w+'.csv')
            
            # use this filepath for river basin load
            # df1.to_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_'+str(v)+'/Avg_Power_amazon/Case_'+c+'/Case'+c+w+'.csv')
            
            # use this filepath for tide load
            # df1.to_csv('/Users/emmadevin/School/ViscousDissipation/Data_Rescaled/Viscosity_'+str(v)+'/Avg_Power_tides/Case_'+c+'/Case'+c+w+'.csv')
        
        v+=1