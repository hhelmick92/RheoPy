# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 08:28:10 2023

@author: hhelmick
"""

import pandas as pd
import numpy as np
from numpy import cos as cos
from numpy import sin as sin
import matplotlib.pyplot as plt
plt.style.use('seaborn')
import seaborn as sns
from scipy.optimize import curve_fit
from scipy import stats
from scipy import optimize
from scipy.interpolate import splrep, splev
from scipy.stats import linregress

# In[]

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

# In[]

def get_amp(df_in, Graph = False, spline_smoothing = 10):
    '''

    Parameters
    ----------
    df_in : Pandas Data Frame
        Data that you want to calculate the roots of
    graph : Boolean
        Do you want to plot your data or not

    Returns
    -------
    list
        [0] amplitude of the first peak
        [1] amplitude of the second peak
        [2] frequency of wave in Hz
        [3] frequency of wave in rad /s 
        [4] the number of roots that are found

    '''
    
    x = df_in['Step time'].astype('float') # read in time data and convert it to a float
    y = df_in['Strain'].astype('float') # read in strian data and convert it to a float
            
    # make a data frame, take the derivative
    df1 = pd.DataFrame() 
    df1['time'] = x
    df1['strain'] = y
    df1['time_diff'] = df1['time'].diff()
    df1['strain_diff'] = df1['strain'].diff()
    df1['deriv'] = df1['strain_diff'] / df1['time_diff']
    
    result = linregress(df1.index, y)
    df1['strain_diff_corr'] = df1['strain_diff'] - result.slope
    df1['deriv_corr'] = df1['strain_diff_corr'] / df1['time_diff']

    # spline the data to get rid of noise that happens from having multiple nearby roots
    bspl = splrep(df1['time'][1:], df1['deriv_corr'][1:], s = spline_smoothing)
    df1['spline'] = splev(x, bspl)
    
    # ignore the first 20 data points, since there is not always a root at T= 0, but if there is , we don't care about it
    cut = 50
    df1 = df1[cut:]
    # get the roots of the splined data
    A = df1['spline'].to_numpy()    
    ch = np.where(A[:-1]*A[1:]<0)[0]
    ch = ch + cut
    
          
    # makes a graph if boolean is true
    if Graph == True:
        # plot the derivative function in blue with the spline over lid in red, 0 grey
        z = np.zeros(len(df1['time']))
        fig, ax = plt.subplots(2,1)
        fig.tight_layout(h_pad = 2)
        ax[0].plot(df1['time'], df1['deriv_corr'], color = 'blue', linewidth = 6.0)
        ax[0].plot(df1['time'], z, color = 'grey')
        ax[0].plot(df1['time'], df1['spline'], color = 'red', linestyle = '-')
        ax[0].set_title('Derivative of strain wrt time - slope')
        ax[0].set_ylabel('del strain / del time')
        
        # plot the second axis showing where the roots line up on the curve
        ax[1].plot(df1['time'], df1['strain'])
        for c in list(range(len(ch))):
            ax[1].plot(df1['time'][ch[c]], df1['strain'][ch[c]], 'ro')
        
        ax[1].set_title('Roots Plotted on Strain Data')
        ax[1].set_xlabel('Time (s)')
        ax[1].set_ylabel('Strain (%)')
        ax[1].plot(df1['time'], result.slope * df1.index + result.intercept)

    else: 
        pass
    
    '''
    Amplitude calculation points taken from
    https://ewoldt.mechanical.illinois.edu/files/2022/10/RHE.Z1-Ewoldt-McKinley-2007-Rheol-Bull-PREPRINT-CreepRinging-FULL-manuscript.pdf
    '''
    
    #calculate the amplitude and frequency based on the roots
    amp1 = df1['strain'][ch[0]] - df1['strain'][ch[1]]
    amp2 = df1['strain'][ch[2]] - df1['strain'][ch[1]]
    freq_hz = 1/ (df1['time'][ch[2]] - df1['time'][ch[0]] )
    rad_s = freq_hz * 2*np.pi
    
    j1 = df1['strain'][ch[0]]
    j2 = df1['strain'][ch[1]]
    j3 = df1['strain'][ch[2]]
    j4 = df1['strain'][ch[3]]
    
    del_new = 2*np.log((j1 - (2*j2) + j3)/(-j2 + (2*j3) - j4))        
    
    return [amp1, amp2, freq_hz, rad_s, len(ch), del_new]

def short_t_fit(xData, a_s):
    '''
    function for fitting a pure quadratic at shor time frames and get a_sm solved
    '''
    strain = (stress / (2*a_s)) * xData**2
    return strain

def Struik(df_in, start_in, end_in, Graph = False, spline_smoothing = 1):
    '''
    Paramters in:
        df_in = data frame that you want to calculate G' and G" on based on the Struik method
        start_in = the starting position of the pure quadratic fit in the fucntion
        end_in = the ending position of where to fit the quadartic function
        
    Parmaeters out:
        [0] the log decrement function delta, calculated based on strain data
        [1] the estimated G' of the material based on the creep rining
        [2] the etimated G" of the material based on the creep ringing data
        [3] the estimated tan delta of the material based on the Struik calculation
        [4] the direct calculation of tan delta as G" / G'
        [5] boolean, checks that the log decrement is < 2 pi, which is a condition of using this method
        [6] the number of roots that were found, it's a sanity check thing
        
    '''
    
    # load in x, y data, make floats
    x = df_in['Step time'].astype('float')
    y = df_in['Strain'].astype('float')
    
    # define where to start and stop the curve fit for finding a_sm
    start = start_in
    end = end_in
    # find a_sm from a curve fit function
    a_sm, error = curve_fit(short_t_fit, x[start:end], y[start:end], bounds = (0, 10))
    
    # plot that curve fit if you want
    if Graph == True:
        short_t_solve = short_t_fit(x[start:end], a_sm[0])
        plt.plot(x[start:end], short_t_solve[start:end])
        plt.plot(x[start:end], y[start:end])
    
    # calculate machine intertia based on the the curve fitted a_sm
    I = b*a_sm
    I = I[0]
    
    # get the roots of the function using the get_amp function
    roots = get_amp(df_in, Graph, spline_smoothing)    
    
    freq = roots[3]
    
    if freq > 50:
        roots = get_amp(df_in, Graph, spline_smoothing = 15)
        
        freq = roots[3]
        del_new = roots[5]
        g_strain = ((I*freq**2) / b) * (1+ ((del_new /(2*np.pi))**2))
        g2_strain = ((I*freq**2)/b) * (del_new/np.pi)
        tan_delta = (del_new / np.pi) * (1+ (del_new/(2*np.pi))**2)
        tan_delta_direct = g2_strain / g_strain
        check = del_new < 2*np.pi
    
    else:
        roots = get_amp(df_in, Graph, spline_smoothing)    
        freq = roots[3]
        del_new = roots[5]
        g_strain = ((I*freq**2) / b) * (1+ ((del_new /(2*np.pi))**2))
        g2_strain = ((I*freq**2)/b) * (del_new/np.pi)
        tan_delta = (del_new / np.pi) * (1+ (del_new/(2*np.pi))**2)
        tan_delta_direct = g2_strain / g_strain
        check = del_new < 2*np.pi
        
    return [del_new, g_strain, g2_strain, tan_delta, tan_delta_direct, a_sm[0], freq, check, roots[4]]

# In[]

def jeff_grabber(df_in, Graph = False, start_in = 0, end_in = 60, omega_in = None):
    '''
    Parameters in:
        df_in = pandas data frame
            this is the data that contains the step teim and strain data that you will fit using the Jeffery Model
        struik_in = list
            This is the data that came from the Struik data, specifically, we need the a_sm value that comes from that fucntion
        true freq = Boolean
            use struik calculated frequncy for True or curve fit it for False
        
    Returns:
        list
        [0] = r2 of the fit for the Jeff model
        [1] = G' estimated from the Jeff Model
        [2] = g" estimated from the Jeff Model
        [3] = frequency from curve fitting
        [4] = check that g > g_critical
        
    Useful publications that I'm following
    https://pubs.rsc.org/en/content/articlehtml/2011/sm/c1sm05399j#cit39    
    https://ewoldt.mechanical.illinois.edu/files/2022/10/RHE.Z1-Ewoldt-McKinley-2007-Rheol-Bull-PREPRINT-CreepRinging-FULL-manuscript.pdf
    
    '''
        
    x = df_in['Step time'].astype('float')
    y = df_in['Strain'].astype('float')
    
    # define where to start and stop the curve fit for finding a_sm
    start = start_in
    end = end_in
    # find a_sm from a curve fit function
    a_sm, error = curve_fit(short_t_fit, x[start:end], y[start:end], bounds = (0, 10))
    
    def Jeff(xData, mus, gs, ns):
        
        A = ( (gs/ns) + (mus / a_sm) ) / (2* (1 + (mus/ns)))
        B = (stress / gs) * ((ns+mus) / ns) * (((2*A*a_sm) / ns) -1)
        if omega_in > 0:
            omega = omega_in
        
        else:
            omega = np.sqrt( ((gs / a_sm) * (ns / (mus + ns))) - (A**2) ) 
        
        
        strain = ((stress / ns)* xData ) - B + np.exp( (-A*xData)) * ((B*cos(omega*xData)) + ((A / omega) * (B - (stress / (ns*A))) * (sin(omega*xData)))) 
        
        return strain
    
    def Jeff_long_return(xData, mus, gs, ns):
        
        A = ( (gs/ns) + (mus / a_sm) ) / (2* (1 + (mus/ns)))
        B = (stress / gs) * ((ns+mus) / ns) * (((2*A*a_sm) / ns) -1)
        B = (stress / gs) * ((ns+mus) / ns) * (((2*A*a_sm) / ns) -1)
        if omega_in > 0:
            omega = omega_in
        
        else:
            omega = np.sqrt( ((gs / a_sm) * (ns / (mus + ns))) - (A**2) ) 
            
        strain = ((stress / ns)* xData ) - B + np.exp( (-A*xData)) * ((B*cos(omega*xData)) + ((A / omega) * (B - (stress / (ns*A))) * (sin(omega*xData))))
        
        t3 = (gs / a_sm) * (ns / (mus + ns))

        return [strain, A, B, omega, t3]

    jeff_params, error = curve_fit(Jeff, x, y, p0 = [1, 10, 10], bounds = ([0,0,0], [10000,10000,10000]))
    jeff_long = Jeff_long_return(x, jeff_params[0], jeff_params[1], jeff_params[2])

    mu = jeff_params[0]
    g = jeff_params[1]
    n = jeff_params[2]

    jeff_solve = Jeff(x, mu, g, n)

    lam1 = (mu+n) / g
    lam2 = n/g

    freq_t1 = g / a_sm
    freq_t2 = n / (mu + n)
    A = ((g / n) + (mu / a_sm)) / (2* (1+(mu/n)))
    freq_solve = np.sqrt(freq_t1 * freq_t2 - A**2)

    jeff_g1 = g * (((lam2 * freq_solve)**2) / (1+ (lam1 * freq_solve)**2))
    jeff_g2 = g * (((lam2*freq_solve) * (1+((lam1**2 - (lam1*lam2))* freq_solve**2))) / (1+ (lam1 * freq_solve)**2))

    g_crit = A**2 * a_sm * (1 + (mu / n))
    g_solved = g

    check = g_solved > g_crit
    
    r2 = jeff_solve.corr(y)**2
    
    if Graph == True:
        plt.plot(x, jeff_solve, linestyle = ':', label = str(round(r2, 3)))
        plt.plot(x, y)
        plt.legend()
    
    return [r2, jeff_g1[0], jeff_g2[0], jeff_long[3], check, lam1, lam2, mu, g, n, A[0]]

# In[]

def BurgerModel(xData_in, G1, N1, G0, N0):
    y = (1/G0) + ((1/G1)*(1-np.exp(-xData_in * G1 / N1))) + (xData_in / N0)
    return y

def BurgerExtract(df_in, color = 'indianred', Graph = False, print_info = False):
    
    x = df_in['Step time'].astype('float')
    y = df_in['Strain'].astype('float')
    y = (y / stress)
    
    burger_params, errors = curve_fit(BurgerModel, x, y)
    
    G1_burger = burger_params[0]
    N1_burger = burger_params[1]
    G0_burger = burger_params[2]
    N0_burger = burger_params[3]
    
    burgerSolved = BurgerModel(x, G1_burger,N1_burger,G0_burger,N0_burger)
    r2 = (y.corr(burgerSolved))**2

    # this is the slope of the end of the creep curve, should compare to N0
    result = stats.linregress(x[450:], y[450:])
    zsv = result[0]**-1

    if Graph == True:
        plt.style.use('seaborn')
        plt.plot(x, burgerSolved, color = color, label = 'fit: R2 ' + str(round(r2,3)), linestyle = 'dashed')
        plt.plot(x, y, color = color, label = 'original data')
        #plt.yscale('log')
        #plt.xscale('log')
        plt.legend(bbox_to_anchor = (1,1))
        plt.xlabel('Time (s)')
        plt.ylabel(('J (1/MPa)'))
    else:
        pass
        
    if print_info == True:
        print("Burger Model r2 value")
        print(r2)
        print('First Elastic component')
        print(G0_burger)
        print('Second Elastic component')
        print(G1_burger)
        print('First Viscous component')
        print(N0_burger)
        print('Second Viscous component')
        print(N1_burger)
    else:
        pass

    return [r2, G0_burger, G1_burger, N0_burger, N1_burger, zsv]

# In[]

def recoverStrain(column_in):
    maxStrain = column_in.max()
    inter1 = column_in[(len(column_in) - 50) : ]
    minStrain = min(inter1)
    
    return ((maxStrain - minStrain) / maxStrain) * 100

def J_MinMax(column_in):
    JMax = column_in.max()
    inter1 = column_in[(len(column_in) - 50) : ]
    JInf = min(inter1)
    
    return [JMax, JInf]

def recoveryModel(df_in, Graph = False, color = 'indianred'):

    x = df_in['Step time'].astype('float')
    strain = df_in['Strain'].astype('float')
    comp = ((strain) / stress) * 1000000
    
    recovered_strain = recoverStrain(strain)
    
    min_max = J_MinMax(comp)
    JMax = min_max[0]
    JInf = min_max[1]
    
    jkv_direct = JMax - JInf
    
    firstEle = 100
    lastEle = 483
    
    x_cut = x[firstEle:lastEle]
    strain_cut = strain[firstEle:lastEle]
    comp_cut = comp[firstEle:lastEle]
    
    def recover(xData_in, JKV_in, B, C):
        y = JInf + (JKV_in *np.exp(-B*xData_in**C))
        return y

    recover_params, _ = curve_fit(recover, x_cut, comp_cut, bounds = ([0, 0,0],[JMax, 10, 10]))

    JKV_solved = recover_params[0]
    B = recover_params[1]
    C = recover_params[2]

    recoverSolved = recover(x, JKV_solved, B,C)
    
    JSM_solved = JMax - (JInf + JKV_solved)
    JSM_no_solve = JMax - (JInf + jkv_direct)

    r2 = (comp.corr(recoverSolved))**2

    if Graph == True:
        plt.plot(x, recoverSolved, linestyle = 'dashed', color = color, label = 'Fit: R2 ' + str(round(r2, 3)))
        plt.plot(x,comp, color = color, label = 'Original data')
        plt.legend(bbox_to_anchor = (1,1))
        plt.xlabel('Time (s)')
        plt.ylabel('J(t)')
    else:
        pass

    return [recovered_strain, JMax, JInf,  jkv_direct, r2, JKV_solved, B, C, JSM_solved, JSM_no_solve]

# In[]

def gibbs_activation(df_in):
    '''
    I'M PRETTY SYRE THE UNITS ON THIS ARE JOULES KELVIN /S'
    '''
    k = 1.380649e-23 # J / K
    h = 6.62607015e-34 # J * s
    T = 298.15 # K
    df_in['1/relax'] = 1/ df_in['relaxation time']

    kTh = (k*T / h)
    kT = k*T
    
    df_in['Gibbs'] = -(np.log((df_in['1/relax'] / kTh) * kT))
    
    arr = df_in['Gibbs'].to_numpy() # J K / s
    
    return arr

def bond_count(df_in):
    
    # assuming probability = 1
    P = 0.01
    df_in['bond count'] = (np.log(df_in['relaxation time']*P)) / np.log(0.5)
    
    arr = df_in['bond count'].to_numpy()
    
    return arr

# In[]

def exp_fit(c1, c2, length = 600):
    
    xData = np.arange(0, length, 0.1)
    def func(x, m, t):
        return m * np.exp(-t * x)

    initialParameters = np.array([1,1])
    fittedParameters, pcov = curve_fit(func, c1, c2, initialParameters)

    expand = func(xData, *fittedParameters)

    return expand
    
def intersect(exp1, exp2):
    
    f = exp1
    g = exp2
    idx = np.argwhere(np.diff(np.sign(f - g))).flatten() 

    return idx

def intersec_plot(x1, y1, y2, length = 600, Graph = False):
    
    y1_exp = exp_fit(x1, y1, length)
    y2_exp = exp_fit(x1, y2, length)
    
    intersection = intersect(y2_exp, y1_exp)
    
    xData = np.arange(0, length, 0.1)
    strain = xData[intersection]
    
    if Graph == True:
        plt.scatter(xData[intersection], y1_exp[intersection], color = 'red', s = 100, zorder = 5, label = 'Cross Over Strain')
        plt.scatter(x1, y1)
        plt.scatter(x1, y2)
        
        plt.plot(xData[2:12], y1_exp[2:12], color = 'red', linewidth = 4, label = 'Linear Region for Mean')
        plt.plot(xData[2:12], y2_exp[2:12], color = 'red', linewidth = 4)
        
        plt.plot(xData, y1_exp)
        plt.plot(xData, y2_exp)
        plt.yscale('log')
        plt.xscale('log')
    
    return [strain, y1_exp, y2_exp, intersection]

def get_cross_over_strain(x1, y1, y2, length = 600):
    
    y1_exp = exp_fit(x1, y1, length)
    y2_exp = exp_fit(x1, y2, length)
    
    intersection = intersect(y2_exp, y1_exp)
    
    xData = np.arange(0, length, 0.1)
    strain = xData[intersection]
        
    return strain





