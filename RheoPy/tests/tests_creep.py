# In[]

# import needed packages

# operating libraries
import glob 

# basic handling libraries
import pandas as pd
import numpy as np

# home grown functions 
from RheoPy import creep

# In[]

# path to creep file
df = pd.read_csv(r'C:\Users\hdhel\Desktop\Coding_Projects\RheoPy\docs\data\pea_xanthan_emulsion.csv', skiprows = 11)

# get the units row. Standard in TRIOS outputs
units = df.iloc[0]

# drop that from the df 
df = df.drop(0)

# convert to numeric, ignore some random errors
df2 = pd.DataFrame()
for c in df:
    df2[c] = pd.to_numeric(df[c], errors = 'coerce')

# get the index where the expriment restarted from creep to compliance
t0 = df2.loc[df2['Step time'] == 0]
idxs = t0.index

# get the two tables for seperate handling
creep_data = df2[0: idxs[1] -1]
recovery = df2[idxs[1] -1 :]

# APPLIED STRESS
stress = 1e-07 #MPa
stress = stress * 1000000 # Pa
######### Need to multiply by 100 for dimensional consistency ########
stress = stress * 100

x = creep_data['Step time'].astype('float')
x = x.dropna(how = 'any')
y = creep_data['Strain'].astype('float')
y = y.dropna(how = 'any')
y = (y / stress)

burger_fit = creep.FitBurgerModel(x, y, Graph = True)
