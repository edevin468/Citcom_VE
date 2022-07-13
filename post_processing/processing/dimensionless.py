#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 10:39:31 2021

Takes integrated dissipation from *.time files and converts it to spatially averaged dissipation in units of mu^2/eta


@author: emmadevin
"""

import numpy as np
import pandas as pd

   
DISS_REAL = []
WORK_REAL = []
ELASTIC_REAL = []
PERIOD =[]


case = '900'
wav = ['', 'A','B','C','D']
wavelengths = {'':2.5, 'A':1.25,'B':5.0,'C':0.625,'D':0.3125}
D =1

for w in wav:
    
    area = wavelengths[w]*D
    volume = wavelengths[w]*D*1
    
    k = 1
    while k < 24:
        
        if len(str(k))== 1:
            i = 'c0'+str(k)
        else:
            i = 'c'+str(k)
            
        
        k+=1
        
        FILE =  '/Users/emmadevin/School/ViscousDissipation/Data_Raw/CASE'+case+w+'/'+i+'.time'
        
        # FILE = '/Users/emmadevin/School/ViscousDissipation/Data//CASE500P1T/'+i+'.time'
        
        df = pd.read_table(FILE, sep="\s+", skiprows=0, header = None)
        df.columns = ['Timestep','Time', 'Total Dissipation','Total Elastic','Total Elastic+Dissipation','Total Work', 'Dissipation', 'Elastic Energy','Work']
    
        
        t = np.array(df['Time'])
        dissipation = df['Dissipation']
        elastic = df['Elastic Energy']
        work = df['Work']
        
        diss_vol = dissipation/(area)
        
        work_vol = work/wavelengths[w]
        
        elastic_vol = elastic/(area)
        
        time = t

        
        df1 = pd.DataFrame()
    
        df1['Time'] = time
        df1['Dissipation'] = diss_vol
        df1['Work'] = work_vol
        df1['Elastic'] = elastic_vol
        
        df1.to_csv('/Users/emmadevin/School/ViscousDissipation/Data_Dimensionless/Total_Dissipation/Case_'+case+'/CASE'+case+w+i+'.csv')
        
        
        
        
        
        
        
        