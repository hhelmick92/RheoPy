
# depdendent libraries of RheoPy
# basic handling libraries
import pandas as pd
import numpy as np

# plotting libraries
import seaborn as sns
import matplotlib.pyplot as plt

# home grown functions 
from RheoPy import creep

# think that this will just be a depdendent
from scipy.optimize import curve_fit
from scipy import stats
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
    y = (1/G0) + ((1/G1)*(1-np.exp(-xData_in * G1 / N1))) + (xData_in / N0)
    return y

def FitBurgerModel(x, y, Graph = False, print_info = False):
        
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

    # this is the slope of the end of the creep curve, should compare to N0
    result = stats.linregress(x[450:], y[450:])
    zsv = result[0]**-1

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
                            name='fit data' + '\n' +
                                  'R2: {}'.format(str(round(r2,3)))
                                ))

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

    return [r2, G0_burger, G1_burger, N0_burger, N1_burger, zsv]

