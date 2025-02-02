# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 08:50:13 2023

@author: hhelmick
"""

import pandas as pd
import glob

from my_functions import get_cross_over_strain
from my_functions import intersec_plot

# In[]

# this is the path in 
path = glob.glob(r'L:\Kokini Lab\Natalia Rodriguez\Spring 2023\rheologicla_data\0_week_saos\1_PH8\*.csv')

# this is the path out
path_out = r'L:\Kokini Lab\Natalia Rodriguez\Spring 2023\rheologicla_data\0_week_saos\SAOS_out\1_PH8.csv'

# In[]

name = []
for fn in path:
    t1 = fn.split('\\')
    t2 = t1[-1]
    t3 = t2[:-4]
    name.append(t3)

# In[]

g1 = []
g2 = []
strain = []

for fn in path:
    df = pd.read_csv(fn, skiprows = 11)
    df = df.drop(0)
    df = df[0:20] ########## BE CARFEUL ON THIS LINE IF YOU CHANGE TIME OF CURVE #########
    for col in df.columns:
        if 'Oscillation strain' in col:
            strain.append(df['Oscillation strain'].astype('float'))
        elif 'Storage modulus' in col:
            g1.append(df['Storage modulus'].astype('float'))
        elif 'time' in col:
            g2.append(df['Loss modulus'].astype('float'))

# In[]

g1_df = pd.DataFrame(g1).transpose()
g1_df.columns = name

g1_linear = g1_df[2:8]
g1_linear = g1_linear * 1000000

g2_df = pd.DataFrame(g2).transpose()
g2_df.columns = name

g2_linear = g2_df[2:8]
g2_linear = g2_linear * 1000000

strain_df = pd.DataFrame(strain).transpose()
strain_df.columns = name

# In[]

g1_linear_means = []
g2_linear_means = []

for c in strain_df.columns:
    g1_linear_means.append(g1_linear[c].mean())
    g2_linear_means.append(g2_linear[c].mean())
    
result_df = pd.DataFrame(g1_linear_means, columns = ['g1_linear_mean'])
result_df['g2_linear_mean'] = g2_linear_means

# In[]

crosses = []
for c in strain_df.columns:
    r = get_cross_over_strain(
        strain_df[c],
        g1_df[c],
        g2_df[c]
        )
    crosses.append(r[0])

result_df['cross_over_strain'] = crosses
result_df.index = name
result_df['linear_tan_delta'] = result_df['g2_linear_mean'] / result_df['g1_linear_mean']

result_df.to_csv(path_out)

# In[PROOF OF CONCEPT PLOT]

import matplotlib.pyplot as plt

cols = strain_df.columns
plot_col = 19
print(cols[plot_col])

concept = intersec_plot(strain_df[cols[plot_col]],
              g1_df[cols[plot_col]],
              g2_df[cols[plot_col]],
              Graph = True
              )

plt.title('Cross over strain figure example')
plt.ylabel("G', G'' (log MPa)")
plt.xlabel('Strain Rate (%)')
plt.legend(loc = 'lower left')

plt.savefig(r'L:\Kokini Lab\Harrison Helmick\Pea Protein\11_emulsions_round2\writing_a_paper\images\SAOS_cross_over_example.png', bbox_inches= 'tight', dpi = 400)


















