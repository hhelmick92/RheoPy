# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 09:32:42 2023

@author: hhelmick
"""

import pandas as pd
import glob
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

# In[READ IN DATA AND CONCATENATE DIFFERENT MODEL OUTPUT FILES]

# READ IN THE RIGHT FILE PATH
path_in = r'L:\Kokini Lab\Natalia Rodriguez\Spring 2023\rheologicla_data\0_week_saos\SAOS_out'
path = glob.glob(path_in + '\*.csv')

# READ EACH FILE AND MAKE A LIST OF DATA FRAMES
l1 = []
for fn in path:
    df = pd.read_csv(fn).transpose()
    l1.append(df.iloc[:,:])

# CONCATENATE THOSE DATA FRAMES INTO A SINGLE DF
df = pd.concat([x for x in l1], axis = 1)

# SOME DATA CLEANUP THINGS, REDINDIXING, MAKING FLOATS, ETC
df2 = df.transpose()
df2.index = df2['Unnamed: 0']
df2 = df2.drop('Unnamed: 0', axis = 1)
df2 = df2.astype('float')

# GENERATING CALLABLES FROM THE INDEX. I.E., GET OUT ALL OF THE TREATMENT, PHI, PH, ETC. DATA
l1 = df2.index

treat = []
conc = []
xg = []
phi = []
ph = []
test = []
rep = []

for ele in l1:
    t1 = ele.split('_')
    treat.append(t1[0])
    conc.append(t1[1])
    xg.append(t1[2])
    phi.append(t1[3])               
    ph.append(t1[4])               
    test.append(t1[5])               
    rep.append(t1[6])               

df2['treat'] = treat
df2['ph'] = ph
df2['rep'] = rep
df2['conc'] = conc
df2['test'] = test

df2['treat_ph'] =  df2['treat'] + '_' + df2['ph'] 
df2['treat_ph'] = df2['treat_ph'].apply(lambda x: x.replace('pH', ''))

df2['conc_ph_treat'] = df2['conc'] + df2['ph'] + '_' + df2['treat']
df2['conc_ph_treat'] = df2['conc_ph_treat'].apply(lambda x: x.replace('proteinpH', '_'))

# READ IN THE CHEMO_RHEO MERGE FILE THATS MADE IN THE CONCAT_CHEM_RHEOLOGICAL SCRIPT
chem_rheo = pd.read_csv(r'L:\Kokini Lab\Natalia Rodriguez\Spring 2023\correlation_sandbox\processed_files\chem_rheo_merge.csv')

chem2 = pd.read_csv(r'L:\Kokini Lab\Natalia Rodriguez\Spring 2023\clean_stat_processed_data\chem_data2.csv')
chem2 = chem2[['treat_ph', 'zp', 'sol', 'sNot', 'phobic_0']]
chem2['treat_ph']

# EXTRACT THE FRESH DATA, SINCE WE ARE ONLY COMPARING FRESH IN THESE CIRCUMSTANCES
fresh = chem_rheo.loc[chem_rheo['time'] == 0]
merge1 = df2.merge(fresh, on = 'conc_ph_treat')
merge = merge1.merge(chem2, on = 'treat_ph')

# CALCULATE THE RMSE OF THE DIFFENT METHODS OF ESTIMAING G1 AND G2
merge['SAOS_ERROR_jeff_g1'] = np.sqrt((merge['g1_linear_mean'] - merge['jeff g1 Pa'])**2)
merge['SAOS_ERROR_jeff_g2'] = np.sqrt((merge['g2_linear_mean'] - merge['jeff g2 Pa'])**2)

merge['SAOS_PERCENT_ERROR_jeff_g1'] = merge['SAOS_ERROR_jeff_g1'] / merge['jeff g1 Pa']
merge['SAOS_PERCENT_ERROR_jeff_g2'] = merge['SAOS_ERROR_jeff_g2'] / merge['jeff g2 Pa']

merge['SAOS_ERROR_Struik_g1'] = np.sqrt((merge['g1_linear_mean'] - merge['struik g1 Pa'])**2)
merge['SAOS_ERROR_Struik_g2'] = np.sqrt((merge['g2_linear_mean'] - merge['struik g2 Pa'])**2)

merge['SAOS_PERCENT_ERROR_Struik_g1'] = merge['SAOS_ERROR_Struik_g1'] / merge['struik g1 Pa']
merge['SAOS_PERCENT_ERROR_Struik_g2'] = merge['SAOS_ERROR_Struik_g2'] / merge['struik g2 Pa']

# In[]

merge = merge[[
# these are the identifiers
        'file_name', 'conc_ph_treat','conc_ph_treat_time', 'conc_x', 'ph_x', 'treat_x', 'treat_ph', 'time', 'test',
# these are the SAOS paremeters           
       'g1_linear_mean', 'g2_linear_mean', 'cross_over_strain', 'linear_tan_delta',
# this is the SAOS error      
       'SAOS_ERROR_jeff_g1', 'SAOS_ERROR_jeff_g2','SAOS_ERROR_Struik_g1', 'SAOS_ERROR_Struik_g2', 
       'SAOS_PERCENT_ERROR_jeff_g1', 'SAOS_PERCENT_ERROR_jeff_g2', 'SAOS_PERCENT_ERROR_Struik_g1', 'SAOS_PERCENT_ERROR_Struik_g2',
# this is the stuff from the struik fitting    
       'log decrement', 'struik g1 Pa', 'struik g2 Pa', 'struik tan del solve', 'struik tan del direct', 'a_sm', 'freq rad per s', 'number of roots', 
# this is from the Jeff model    
       'r2_jeff', 'jeff g1 Pa', 'jeff g2 Pa', 'jeff freq', 'relaxation time', 'retardation time', 'jeff mu', 'jeff g', 'jeff n', 'jeff damping', 
# this is from the Burger model     
       'r2_burger', 'G0', 'G1', 'N0', 'N1', 'burger retardation time', 'no / go', 'go / no', 'g1 / n1', '1/relax',
# this is from the recovery model    
       'r2_recover', 'zsv', 'recoverStrain', 'JMax', 'JInf', 'jkv_direct', 'JKV_solved', 'B', 'C', 'JSM_solved', 'JSM_no_solve',
# this is a list of molecular estimates
       'Gibbs', 'gibbs_act', 'bond count', 'g1_est', 
# this is the physicochemcial data        
       'zp_x', 'zp std', 'zp tuk', 'sol_x', 'sol std', 'sol tuk', 'sNot_y', 'sNot std', 'sNot tuk', 'sNot percent increase', 'abs_zp', 'phobic_0', 
       ]]

# In[]
merge.columns = [[
# these are the identifiers
        'file_name', 'conc_ph_treat','conc_ph_treat_time', 'conc', 'ph', 'treat', 'treat_ph', 'time2', 'test',
# these are the SAOS paremeters           
       'g1_linear_mean', 'g2_linear_mean', 'cross_over_strain', 'linear_tan_delta',
# this is the SAOS error      
       'SAOS_ERROR_jeff_g1', 'SAOS_ERROR_jeff_g2','SAOS_ERROR_Struik_g1', 'SAOS_ERROR_Struik_g2', 
       'SAOS_PERCENT_ERROR_jeff_g1', 'SAOS_PERCENT_ERROR_jeff_g2', 'SAOS_PERCENT_ERROR_Struik_g1', 'SAOS_PERCENT_ERROR_Struik_g2',
# this is the stuff from the struik fitting    
       'log decrement', 'struik g1 Pa', 'struik g2 Pa', 'struik tan del solve', 'struik tan del direct', 'a_sm', 'freq rad per s', 'number of roots', 
# this is from the Jeff model    
       'r2_jeff', 'jeff g1 Pa', 'jeff g2 Pa', 'jeff freq', 'relaxation time', 'retardation time', 'jeff mu', 'jeff g', 'jeff n', 'jeff damping', 
# this is from the Burger model     
       'r2_burger', 'G0', 'G1', 'N0', 'N1', 'burger retardation time', 'no / go', 'go / no', 'g1 / n1', '1/relax',
# this is from the recovery model    
       'r2_recover', 'zsv', 'recoverStrain', 'JMax', 'JInf', 'jkv_direct', 'JKV_solved', 'B', 'C', 'JSM_solved', 'JSM_no_solve',
# this is a list of molecular estimates
       'Gibbs', 'gibbs_act', 'bond count', 'g1_est', 
# this is the physicochemcial data        
       'zp', 'zp std', 'zp tuk', 'sol', 'sol std', 'sol tuk', 'sNot', 'sNot std', 'sNot tuk', 'sNot percent increase', 'abs_zp', 'phobic_0', 
       ]]

merge.to_csv(r'L:\Kokini Lab\Natalia Rodriguez\Spring 2023\correlation_sandbox\processed_files\saos_merged_data.csv')

