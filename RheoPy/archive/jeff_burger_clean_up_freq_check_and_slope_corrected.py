# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 09:17:45 2023

@author: hhelmick
"""

# In[IMPOR THE WORLD]
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import matplotlib.pyplot as plt
plt.style.use('seaborn')

import glob
import time as T

from my_functions import Struik
from my_functions import jeff_grabber
from my_functions import BurgerExtract
from my_functions import recoveryModel

time1 = T.time()

# In[THIS IS THE SECTION THAT YOU CAHNGE TO RUN THIS CODE]

##################################################################################
########## THIS IS THE ONLY PART OF THIS CODE YOU SHOULD NEED TO CHANGE ##########
##################################################################################

# PUT THE PATH TO WHERE YOU WANT TO RRUN THE CODE HERE. THE FILE PATH CANNOT HAVE A \\ AT THE END
path_in = r'/Volumes/ag_fdsc/Labs/Kokini Lab/Natalia Rodriguez/Spring 2023/rheologicla_data/4_week_creep/5_PH3'

# PUT THE PATH WHERE YOU WANT THE RESULTS HERE. THE FILE PATH CANNOT HAVE A \\ AT THE END
path_out = r'L:\Kokini Lab\Natalia Rodriguez\Spring 2023\rheologicla_data\creep_out\model_results_2'

# PUT WHAT YOU WANT TO NAME YOUR FILE HERE. iT SHOULD CONTAIN NO \\ IN THE TREATMENT NAME
treatment = '05_PH3_4week'

##################################################################################
########################## STOP CHANGING THINGS NOW ##############################
##################################################################################

print('##################################################################################')
print('Here we go!')
print('##################################################################################')


# In[READ IN THE FILE PATHS]

# Generate the list of files to read in using Glob
path = glob.glob(path_in + '\*.csv')

# Generate a list of all of the final names, ignoring the rest of the path. These will be the columns in you final data frame
name = []
for fn in path:
    t1 = fn.split('\\')
    t2 = t1[-1]
    t3 = t2[:-4]
    name.append(t3)
    
# before running the code, double check tha your file path is correct.
assert len(name) > 0, 'Your file path is empty! Double check your file in!'

# In[RUN THE STRUIK CALCULATIONS ON THE CREEP RINING DATA]

struik_data = []

for fn in path:
    df = pd.read_csv(fn, skiprows = 11) # skips the first 11 rows that contain meta data
    df = df.drop(0) # drop the units row
    df = df[0:205] # define how many of the data points you will use for the Struik model. This should match the Jeffery part. 
    
    struik_extract = Struik(df, start_in = 0, end_in = 50, Graph = False, spline_smoothing = 1)
    struik_data.append(struik_extract)

struik_df = pd.DataFrame(struik_data).transpose()
struik_df.columns = name
struik_df.index = (['log decrement', 'struik g1 Pa', 'struik g2 Pa', 'struik tan del solve', 'struik tan del direct', 'a_sm', 'freq rad per s', 'log del less than 2 pi', 'number of roots'])

# In[]

jeff_data = []

freq_list = struik_df.loc['freq rad per s']

for fn in zip(path, freq_list):
    df = pd.read_csv(fn[0], skiprows = 11)
    df = df.drop(0)
    df = df[0:205] ########## BE CARFEUL ON THIS LINE, YOU'LL FIT TO WHATEVER THIS LENGTH IS ~ 2-3 SECONDS RIGHT NOW #########
    
    jeff_extract = jeff_grabber(df, Graph = False, start_in = 0, end_in = 50, omega_in = fn[1])
    jeff_data.append(jeff_extract)

jeff_df = pd.DataFrame(jeff_data).transpose()
jeff_df.columns = name
jeff_df.index = (['r2_jeff', 'jeff g1 Pa', 'jeff g2 Pa', 'jeff freq', 'g > g critical','relaxation time', 'retardation time', 'jeff mu', 'jeff g', 'jeff n', 'jeff damping'])

# In[]

burger_data = []    

for fn in path:
    df = pd.read_csv(fn, skiprows = 11)
    df = df.drop(0)
    df = df[0:483] ########## BE CARFEUL ON THIS LINE, YOU'LL FIT TO WHATEVER THIS LENGTH IS ~ 2-3 SECONDS RIGHT NOW #########
    
    burger_temp = BurgerExtract(df, Graph = False)
    burger_data.append(burger_temp)

burger_df = pd.DataFrame(burger_data).transpose()
burger_df.columns = name
burger_df.index = (['r2_burger', 'G0', 'G1', 'N0', 'N1', 'zsv'])

# In[GET RECOVER DATA]
           
recover_data = []    

for fn in path:
    df = pd.read_csv(fn, skiprows = 499)
    df = df.drop(0)
    df = df[0:483] ########## BE CARFEUL ON THIS LINE, YOU'LL FIT TO WHATEVER THIS LENGTH IS ~ 2-3 SECONDS RIGHT NOW #########
    
    recover_temp = recoveryModel(df, Graph = False)
    recover_data.append(recover_temp)

recover_df = pd.DataFrame(recover_data).transpose()
recover_df.columns = name
recover_df.index = (['recoverStrain', 'JMax', 'JInf',  'jkv_direct', 'r2_recover', 'JKV_solved', 'B', 'C', 'JSM_solved', 'JSM_no_solve'])

# In[]

df_out = pd.concat([struik_df, jeff_df, burger_df, recover_df])

df_out.to_csv(path_out + "\\" + treatment + ".csv")

elapsed = round(T.time() - time1, 2)
print('##################################################################################')
print('Holy smokes! It only took {} seconds to process all {} of those creep curves!'.format(elapsed, len(name)))
print('##################################################################################')

































