
# depdendent libraries of RheoPy
# basic handling libraries
import pandas as pd
import numpy as np

# plotting libraries
import seaborn as sns
import matplotlib.pyplot as plt

# think that this will just be a depdendent
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.interpolate import splrep, splev
from scipy.stats import linregress

# plotly dependent
import plotly.express as px
from plotly.offline import plot as show_plotly 
import plotly.graph_objects as go

def test_function(test_string):
    '''
    Function that shows the program is working
    
    Type something and it returns it as a string
    '''
    
    return (str(test_string))

def BurgerModel(xData_in, G1, N1, G0, N0):
    '''
    Equation of the Burger Model. 
    Returns y given X data from a curve fit.
    '''
    y = (1/G0) + ((1/G1)*(1-np.exp(-xData_in * G1 / N1))) + (xData_in / N0)
    return y

def FitBurgerModel(x, y, Graph = False, print_info = False):
    """
    Fits the burger model to the creep data entered into the function.
    Inputs:
        Mandatory:
        x- Time data in seconds from a Creep Experiment
        y- Compliance data in MPa from a creep experiment*

        Optional:
        Graph- Interactive Plotly plot of the data
        print_info- quick view of the curve fit parameters returned
        
    Returns:
        list- r2, G0, G1, N0, N1
        These are the r2 of the curve fit along with the solved paramters from the model. 

    * Compliance = Strain (%) / Stress where stress is a constnat in the expriment. 
    These can be take form the equipment or calculated based on your design. 
    """    
    
    def BurgerModel(xData_in, G1, N1, G0, N0):
        y = (1/G0) + ((1/G1)*(1-np.exp(-xData_in * G1 / N1))) + (xData_in / N0)
        return y
    
    burger_params, errors = curve_fit(BurgerModel, x, y)
    
    G1_burger = burger_params[0]
    N1_burger = burger_params[1]
    G0_burger = burger_params[2]
    N0_burger = burger_params[3]
    
    burgerSolved = BurgerModel(x, G1_burger,N1_burger,G0_burger,N0_burger)
    r2 = (y.corr(burgerSolved))**2

    data = pd.DataFrame()
    data['Step time'] = x
    data['Compliance'] = y
    data['fit'] = burgerSolved
    if Graph == True:

        fig = go.Figure()

        # Add traces
        fig.add_trace(go.Scatter(
                            x=data['Step time'], 
                            y=data['Compliance'],
                            name='original data',
                            mode = 'markers',
                            marker = {'color' : 'black'}
                            ))
        fig.add_trace(go.Scatter(
                            x=data['Step time'], 
                            y=data['fit'],
                            mode='lines',
                            line = {'color' : 'indianred'},
                            name='Burger Model' + '\n' +
                                  'R2: {}'.format(str(round(r2,3)))
                                ))

        fig.update_layout(
            xaxis_title= 'Tims (s)', 
            yaxis_title= 'Compliance (MPa)'
                        )

        show_plotly(fig)

        '''
        # seaborn plot of the data
        plt.plot(x, burgerSolved, color = color, label = 'fit: R2 ' + str(round(r2,3)), linestyle = 'dashed')
        plt.plot(x, y, color = color, label = 'original data')
        #plt.yscale('log')
        #plt.xscale('log')
        plt.legend(bbox_to_anchor = (1,1))
        plt.xlabel('Time (s)')
        plt.ylabel(('J (1/MPa)'))
        '''
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

    return [r2, G0_burger, G1_burger, N0_burger, N1_burger]


# This was in the burger fit. We need a better way to get the slop of the end of the expriemnt
# maybe infer the time step then work back N number of points
    # this is the slope of the end of the creep curve, should compare to N0
    #result = stats.linregress(x[450:], y[450:])
    #zsv = result[0]**-1

def DampingWaveData(x, y, 
             cut_off = 50,
             region_start = 0, region_end = 205, 
             Graph = False, spline_smoothing = 10):
    """
    Extracts relevant features from the ringing region of creep
    Relies on the residual for the derivative to for finding roots
    there are the peaks in thw ave

    Inputs:
        Mandatory:
        x- Time data in seconds from a Creep Experiment
        y- Strain data from a Creep experiment 

        Optional:
        cut_off*- the number of data points to drop at the beginning of the expierment
        region_start- the beginning of your wave. Use the plotly plot from FitBurgerModel to help define this
        region_end- the point where your wave ends
        Graph- plot the derivative and the roots
        spline_smoothing- how tight of a spline to fit to the derivative
        
    Returns:
        list
        [0] Amplitude of the first iteration of the wave
        [1] Amplitude of the second iteration of the wave
        [2] Wave frequency (Hz)
        [3] Wave frequency (rad /s)
        [4] Count of roots in the derivative
        [5] the log decrement function delta

    *every wave will have a root at t=0. So we drop a few of the first data points, since this isn't
    the actual wave
    """    
    
    # cut the data to user defined range of wave
    x = x[region_start:region_end]
    y = y[region_start:region_end]
            
    # make a data frame, take the derivative
    df1 = pd.DataFrame() 
    df1['time'] = x
    df1['strain'] = y
    df1['time_diff'] = df1['time'].diff()
    df1['strain_diff'] = df1['strain'].diff()
    df1['deriv'] = df1['strain_diff'] / df1['time_diff']
    
    # get the residual for plotting getting the roots of the wave
    result = linregress(df1.index, y)
    df1['strain_diff_residual'] = df1['strain_diff'] - result.slope
    df1['deriv_residual'] = df1['strain_diff_residual'] / df1['time_diff']

    # spline the data to get rid of noise that happens from having multiple nearby roots
    bspl = splrep(df1['time'][1:], df1['deriv_residual'][1:], s = spline_smoothing)
    df1['spline'] = splev(x, bspl)
    
    # ignore the first N data points, since there is not always a root at T= 0, but if there is , we don't care about it
    cut = cut_off
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
        ax[0].plot(df1['time'], df1['deriv_residual'], 
                   color = 'blue', linewidth = 6.0)
        ax[0].plot(df1['time'], z, color = 'grey')
        ax[0].plot(df1['time'], df1['spline'], 
                   color = 'red', linestyle = '-')
        ax[0].set_title('Derivative of strain wrt time - slope')
        ax[0].set_ylabel('del strain / del time')
        
        # plot the second axis showing where the roots line up on the curve
        ax[1].plot(df1['time'], df1['strain'])
        for c in list(range(len(ch))):
            ax[1].plot(df1['time'][ch[c]], df1['strain'][ch[c]], 'ro')
        
        ax[1].set_title('Roots Plotted on Strain Data')
        ax[1].set_xlabel('Time (s)')
        ax[1].set_ylabel('Strain (%)')
        ax[1].plot(df1['time'], result.slope * df1.index + result.intercept, label = 'residual')

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

def short_t_fit(xData, a_s, stress):
    '''
    function for fitting a pure quadratic at short time frames and get a_sm solved
    '''
    strain = (stress / (2*a_s)) * xData**2
    return strain

def FitASM(x, y, start_in, end_in, stress, 
           Graph = False):
    '''
    Extracts the asm value using a pure quadratic from the first N points of data
    This paramter is required fro Struik and Jeffrey models

    Inputs:
        Mandatory:
        x- Time data in seconds from a Creep Experiment
        y- Strain data from a Creep experiment 
        start_in- the point where this data shuold start fitting
        end_in- the point where the data should stop fitting
        stress- the constant stress applied on the rheometer
        
    Returns:
        list
        [0] asm estimated from the quadratic curve fit

    * can use the plotly interactive plot from FitBurgerModel to help find where this region is
    '''
    
    # load in x, y data, make floats
    x = x
    y = y
    
    # define where to start and stop the curve fit for finding a_sm
    start = start_in
    end = end_in
    # find a_sm from a curve fit function
    a_sm, error = curve_fit(short_t_fit, x[start:end], y[start:end], bounds = (0, 10))
    
    # plot that curve fit if you want
    if Graph == True:
        short_t_solve = short_t_fit(x[start:end], a_sm[0], stress = stress)
        plt.plot(x[start:end], short_t_solve[start:end])
        plt.plot(x[start:end], y[start:end])

    return a_sm

def Struik(x, y, asm, b, freq, del_new):
    '''    
    Estimates linear G' and G" based on the wave data in the creep ringing region

    Inputs:
        Mandatory:
        x- Time data in seconds from a Creep Experiment
        y- Strain data from a Creep experiment 
        asm- parameter pulled from FitASM
        b- geometric constant from your rheometer
        freq- frequency of the wave in the ringing region. Use DampingWaveData to get this parameter
        del_new- the log decrement function delta. Use DampingWaveData to get this parameter

    Returns:
        list
        [0] the estimated G' of the material based on the creep rining
        [1] the etimated G" of the material based on the creep ringing data
        [2] the estimated tan delta of the material based on the Struik calculation
        [3] the direct calculation of tan delta as G" / G'
        [4] boolean, checks that the log decrement is < 2 pi, which is a condition of using this method
    
    '''
    
    # load in x, y data, make floats
    x = x
    y = y
    
    # calculate machine intertia based on the the curve fitted a_sm
    I = b*asm
    I = I[0]
        
    if freq > 50:
        
        freq = freq
        del_new = del_new
        g_strain = ((I*freq**2) / b) * (1+ ((del_new /(2*np.pi))**2))
        g2_strain = ((I*freq**2)/b) * (del_new/np.pi)
        tan_delta = (del_new / np.pi) * (1+ (del_new/(2*np.pi))**2)
        tan_delta_direct = g2_strain / g_strain
        check = del_new < 2*np.pi
    
    else:
        freq = freq
        del_new = del_new
        g_strain = ((I*freq**2) / b) * (1+ ((del_new /(2*np.pi))**2))
        g2_strain = ((I*freq**2)/b) * (del_new/np.pi)
        tan_delta = (del_new / np.pi) * (1+ (del_new/(2*np.pi))**2)
        tan_delta_direct = g2_strain / g_strain
        check = del_new < 2*np.pi
        
    return [g_strain, g2_strain, tan_delta, tan_delta_direct, check]
