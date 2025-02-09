# In[]

# import needed packages

# operating libraries
import glob 

# basic handling libraries
import pandas as pd
import numpy as np

# home grown functions 
from RheoPy import creep

'''
Example 1- Fit Burger model to creep data
'''

# path to creep file
df = pd.read_csv(r'C:\Users\hdhel\Desktop\Coding_Projects\RheoPy\docs\data\pea_xanthan_emulsion.csv', skiprows = 11)

# get the units row. Standard in TRIOS outputs
units = df.iloc[0]

# drop that from the df 
df = df.drop(0)

# convert to numeric, handle some load errors
df2 = pd.DataFrame()
for c in df:
    df2[c] = pd.to_numeric(df[c], errors = 'coerce')

# get the index where the expriment restarted from creep to compliance
t0 = df2.loc[df2['Step time'] == 0]
idxs = t0.index

# get the two tables for seperate handling
creep_data = df2[0: idxs[1] -1]
recovery = df2[idxs[1] -1 :]

# Enter Rheological constants
# APPLIED STRESS
stress = 1e-07 #MPa
stress = stress * 1000000 # Pa
######### Needed to multiply by 100 for dimensional consistency (%) ########
stress = stress * 100

# clean up the time and compliance data for the model
x = creep_data['Step time'].astype('float')
x = x.dropna(how = 'any')
y = creep_data['Strain'].astype('float')
y = y.dropna(how = 'any')
y = (y / stress)

# get the burger model data from your creep data
burger_fit = creep.FitBurgerModel(x, y, Graph = True)

# this plot is interactice so that you can find the cutoffs
# for the next section

# In[]

# import needed packages

# operating libraries
import glob 

# basic handling libraries
import pandas as pd
import numpy as np

# home grown functions 
from RheoPy import creep

'''
Example 2- Get G'  and G" estimates from Struik's Method
'''

##############################################
# INPUT CONSTANTS FROM YOUR RHEOMETER
##############################################

# APPLIED STRESS
stress = 1e-07 #MPa
stress = stress * 1000000 # Pa
######### FOR SOME REASON THE 1 / MPA VALUE IN TRIOS ONLY MATHCES WHEN MULTIPLYING STRESS BY 100 ########
stress = stress * 100

# GEOMETRY DIMENSIONS
R_mm = 20 # milimeters
R_M = R_mm / 1000 # meters

# GAP HEIGHT
h_micron = 600 # micron
h_M = h_micron / 1000000 # meters

# CALCUATION OF b GEOMETRIC FACTOR FOR A PARALELL PLATE
b = (np.pi * R_M **4) / 2*h_M

##############################################
##############################################

# path to creep file
df = pd.read_csv(r'C:\Users\hdhel\Desktop\Coding_Projects\RheoPy\docs\data\pea_xanthan_emulsion.csv', skiprows = 11)

# get the units row. Standard in TRIOS outputs
units = df.iloc[0]

# drop that from the df 
df = df.drop(0)

# convert to numeric, handle some load errors
df2 = pd.DataFrame()
for c in df:
    df2[c] = pd.to_numeric(df[c], errors = 'coerce')

# get the index where the expriment restarted from creep to compliance
t0 = df2.loc[df2['Step time'] == 0]
idxs = t0.index

# get the two tables for seperate handling
creep_data = df2[0: idxs[1] -1]
recovery = df2[idxs[1] -1 :]

# clean up time, strain data for model 
x = creep_data['Step time'].astype('float')
x = x.dropna(how = 'any')
y = creep_data['Strain'].astype('float')
y = y.dropna(how = 'any')

# get the wave data needed for struik fits
# this function is left seperate because it is
# not specific to creep rining and can be used
# on other dampening wave functions 
wave_data = creep.DampingWaveData(x, y, 
                            spline_smoothing = 1,
                            Graph = True
                        )

# get the required inputs for the struik data
freq = wave_data[3]   
del_new = wave_data[-1]

# this function fits a polynomial to the 
# first N data points. Used to help estimate the
# intertia of the instrument
asm = creep.FitASM(x = x, 
                    y = y, 
                    stress = stress, 
                    start_in = 0,
                    end_in = 50,
                    Graph = False)

# This uses the calcualtions outlined in Struik to
# estimate G', G" from the data
struik_data = creep.Struik(x = x,
                           y = y, 
                           asm = asm,
                           b= b,
                           freq = freq,
                           del_new = del_new,
                           )

# In[]

'''
A few next steps
clean up the graphs to show what the lines are
the orange line in the Roots plotted on starin data is the curve fit 
Fit used to calculate the residual 

The top plot is the residual and the spline fit for it
the blue line is the actual data 
'''

# In[]

