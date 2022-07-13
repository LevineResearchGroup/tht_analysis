#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 01:39:20 2022

@author: bab226
"""

'''Color map'''
import matplotlib.pyplot as plt
import numpy as np
import lmfit
import seaborn as sns
from scipy.stats import linregress
from uncertainties import ufloat
import tht_fitting as fitme
import numdifftools as nd

def get_cmap(n, name='viridis'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

'''Fitting lines'''
def linfit(x, y):
    
    '''Fit for linear regression'''
    slope, intercept, rvalue, pvalue, stderr = linregress(x,y)
    return slope, intercept

'''Calculating intercept'''
def find_nearest(array, value):
    '''Find nearest value in array and return'''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

'''Guess kinetic paramaters'''
def guess_fit(x, y, label, interval, savetitle, headers, i):
    '''
    Parameters
    ----------
    x : df
        time.
    y : df
        intensity.
    label : string
        selection of data.

    Returns
    -------
    tlag, t50, label
    '''
    
    ''' Define cutoffs '''
    cutoff_lag = 0.15
    cutoff_exp_start = 0.3
    cutoff_exp_end = 0.8 - cutoff_exp_start
    
    '''Guess fitting paramaters based on data'''
    r1_guess = np.mean(y[:5])
    r2_guess = np.mean(y[len(y)-5:len(y)])
    ind = fitme.find_nearest(y, r2_guess/2)
    t50_guess = x[ind]
    
    # m1_guess = (y[5]-y[0])/(x[5]-x[0])
    # m2_guess = (y[len(x)-1]-y[len(x)-5])/(x[len(x)-1]-x[len(x)-5])
    
    if (r1_guess > r2_guess) or (r2_guess < 2000):
        print ("Signal is not WORTHY! Check data.")
    
        return np.nan, np.nan, np.nan, np.nan
    
    else:
        
        tmin = int(x[0])
        tmax = int(x[len(x)-1])
        '''Lag phase'''
        y_cut_lag = fitme.find_nearest(y, r2_guess*cutoff_lag)
            
        tmin = int(x[0])
        tmax = int(x[len(x)-1])
        y_lag = y[tmin:y_cut_lag]
        x_lag = x[tmin:y_cut_lag]
        
        m_lag, b_lag = fitme.linfit(x_lag, y_lag)
        
       
        #plt.figure()
        '''Exp phase'''
        t50 = int(x[ind])
        tind = int(t50/interval)
        
        y_exp = y[tind-5:tind+5]
        
        x_exp = x[tind-5:tind+5]
        
        m_exp, b_exp = fitme.linfit(x_exp, y_exp)
    
        '''Equations of lines'''
        y_fit_lag = m_lag*x+b_lag
        y_fit_exp = m_exp*x+b_exp
        
        '''Get intercepts'''
        xi = (b_lag-b_exp) / (m_exp-m_lag)
        yi = m_lag * xi + b_lag
        
        print('tlag %s = %s' %(label, xi))
        tlag_guess = xi
        
        '''Tau guess'''
        tau_guess = (t50_guess - tlag_guess) / 2
        
        '''Plot guess fit'''
        '''Color MAP'''
        color = sns.color_palette('deep', len(headers), desat=1)
        color_ind = list(range(0,len(color), 3))
        step = int(i/3)
        color = color[color_ind[step]]
        
        plt.plot(x, y, label=label, color=color, zorder=1)
        #plt.plot(x_lag, y_lag, x_lag[tmin], color=colors)
        plt.plot(x, y_fit_lag, label="fit1", color=color, linestyle='--', linewidth=0.7, zorder=4)
        #plt.plot(x_exp, y_exp, color=colors)
        plt.plot(x, y_fit_exp, label="fit2", color=color, linestyle='--', linewidth=0.7, zorder=4)
        #plt.axvline(x = xi, color = color, linestyle='--', label = 'tlag', zorder = 5)
        plt.scatter(xi, yi, color='black')
        plt.scatter(t50_guess, y[ind], color='purple')
        
        #plt.legend(loc="best", bbox_to_anchor=(1.3,1.1), fontsize=12)
        plt.xlabel('time (min)', fontsize=16)
        plt.ylabel('Intensity', fontsize=16)
        plt.xlim(tmin,tmax)
        plt.ylim(0,10000)
        
        c = color
        
        '''Max velocity of monomer consumption'''
        vt50 = m_exp
        
        return t50_guess, tau_guess, vt50, c

'''From Dr. Mirnaker's 2002 paper on fitting ThT curves.'''
def miranker_sigmoid(t, m1, m2, r1, r2, t50, tau):
    return(((m2*t + r2) * (1 + np.e**(t50 - t)/tau)**-1) + ((m1*t + r1) * (1 - (1 + np.e**(t50 - t)/tau)**-1)))

'''Least-squares fitting'''
def sigfit(x, y, label, savetitle, i, tau_guess, c):
    '''
    Parameters
    ----------
    x : df
        time.
    y : df
        intensity.
    label : string
        selection of data.
    savetitle : string
        description of sample
    i : int
        keep track of column and plotting colors
    tau_guess : float
        tau guess from guess_fit funciton
    c : color code output from guess_fit (for consistency)
    
    Returns
    -------
    Paramaters from best_fit
    params include:
        r1 + m1*t 
    '''
        
    '''Guess fitting paramaters based on data'''
    r1_guess = np.mean(y[:5])

    m1_guess = np.abs((y[10]-y[0])/(x[10]-x[0]))
    m2_guess = np.abs((y[len(x)-1]-y[len(x)-10])/(x[len(x)-1]-x[len(x)-10]))
    
    if m1_guess == 0:
        m1_guess = 1
    
    if m2_guess == 0:
        m1_guess = 1
        
    max_val = y[len(x)-1]
    tmin = x[0]
    tmax = x[len(x)-1]
    
    r2_guess = np.mean(y[len(y)-5:len(y)]-1)-(tmax*m2_guess)
    
    y50 = abs(y - r2_guess/2)
    ind = np.argmin(y50)
    t50_guess = x[ind]
    #tau_guess = 0.01
    if (r1_guess > r2_guess) or (max_val < 500):
        print ("\n Signal is not WORTHY! Check data. \n")
    
        return np.nan
    
    else:
        '''Best Fitting Method using Least-squares'''
        model = lmfit.Model(miranker_sigmoid)
        params = model.make_params()
        params['r1'].set(min=0, max=max_val, value=r1_guess, vary=True)
        params['r2'].set(min=0, max=max_val, value=r2_guess, vary = True)
        params['t50'].set(min=0, max=tmax, value=t50_guess, vary=True)
        params['m1'].set(min=0, max=m1_guess*100, value=m1_guess, vary=True)
        params['m2'].set(min=0, max=m2_guess*100, value=m2_guess, vary=False)
        params['tau'].set(min=0, max=tau_guess*10, value=tau_guess, vary=True)
        
        print (params)
        #weights = 1/err
        
        fitres = model.fit(y, t=x, params=params, method='least-squares')
        # fitres = model.fit(y, t, params=params, method='least_squares',
        #                    weights=weights[t<t_max])
        print('\nList of fitted parameters for %s: \n' % model.name)
        fitres.params.pretty_print(colwidth=10, columns=['value', 'stderr', 'min', 'max'])
        
        print(fitres.fit_report())
        
        t50 = fitres.values['t50']
        tau = fitres.values['tau']
        tlag = t50 - 2*tau
        
        print('tlag %s = %s' %(label, tlag))
        
        plt.plot(x, fitres.best_fit, label=label, color='black', linestyle='--')
        # plt.legend(loc="best", bbox_to_anchor=(1.2,1.1), fontsize=12)
        plt.xlabel('time (min)')
        plt.ylabel('Intensity')
        plt.xlim(tmin,tmax)
        plt.ylim(0,10000)
        
        return fitres
    
    
    ##################
'''Note: This function is depreciated and replaced by guess_fit and sigfit in tht_analysis_v6.py'''
    ##################

'''Fitting ThT curve'''
def fit_tht(x, y, label, cutoff_lag, cutoff_exp, color, i, savetitle):
    '''

    Parameters
    ----------
    x : df
        time.
    y : df
        intensity.
    label : string
        selection of data.
    cutoff : float
        cutoff for lag and exp selection ranges.

    Returns
    -------
    tlag, t50, kapp, max_int.

    '''
    
    yf_guess = np.mean(y[-5:])
    yf_guess_err = np.std(y[-5:])
    yf_guess = ufloat(yf_guess, yf_guess_err)
    yi_guess = np.mean(y[:5])
    yi_guess_err = np.std(y[:5])
    yi_guess = ufloat(yi_guess, yi_guess_err)

    if (yi_guess.nominal_value > yf_guess.nominal_value) or (yf_guess < 2000):
        print ("Signal decreases OR is not WORTHY! Check data.")
    
        return np.nan, np.nan, np.nan, np.nan, np.nan
    else:
        
        '''Selection cutoff for linear regime of lag phase'''    
        y_cut = find_nearest(y, yf_guess.nominal_value*cutoff_lag)
    
        tmin = int(x[0])
        tmax = int(x[len(x)-1])
        y_lag = y[tmin:y_cut]
        x_lag = x[tmin:y_cut]
      
        m_lag, b_lag, r_lag, p_lag, err_lag = linfit(x_lag, y_lag)
        m_lag = ufloat(m_lag, err_lag)
        y_fit_lag = m_lag.nominal_value*x+b_lag
        
        '''Selection cutoff for linear regime of exponential phase'''
        cutoff_min = cutoff_exp
        cutoff_max = 1 - cutoff_min
    
        y_cut_min = find_nearest(y, yf_guess.nominal_value*cutoff_min)
        y_cut_max = find_nearest(y, yf_guess.nominal_value*cutoff_max)

        y_exp = y[y_cut_min:y_cut_max]
        x_exp = x[y_cut_min:y_cut_max]
           
        m_exp, b_exp, r_exp, p_exp, err_exp = linfit(x_exp, y_exp)
        y_fit_exp = m_exp*x+b_exp
        
        '''Get intercept'''
        xi = (b_lag-b_exp) / (m_exp-m_lag)
        yi = m_lag * xi + b_lag
        print('tlag %s = %s' %(label, xi))
        tlag = xi
        #tlag = int(np.round(xi.nominal_value))
        
        '''Get T50'''
        y_50 = yf_guess/2
        ind_50 = find_nearest(y, y_50)
        t_50 = ufloat(x[ind_50], y_50.std_dev)
        
        '''kapp'''
        tau = (t_50 - tlag)/2
        kapp = 1.0/tau
        
        '''max_int'''
        max_int = yf_guess
        
        '''Plotting data and fits'''
        color_ind = list(range(0,len(color), 3))
        step = int(i/3)
        colors = color[color_ind[step]]
        
        plt.plot(x, y, label=label, color=colors, zorder=1)
        plt.plot(x_lag, y_lag, x_lag[tmin], color=colors)
        plt.plot(x, y_fit_lag, label="fit1", color=colors, linestyle='-', linewidth=0.7, zorder=4)
        plt.plot(x_exp, y_exp, color=colors)
        plt.plot(x, y_fit_exp, label="fit2", color=colors, linestyle='-', linewidth=0.7, zorder=4)
        #plt.axvline(x = xi.nominal_value, color = colors, linestyle='--', label = 'tlag', zorder = 5)
        plt.scatter(xi.nominal_value, yi.nominal_value, color='black')
        plt.scatter(t_50.nominal_value, y_50.nominal_value, color='purple')
        
        plt.legend(loc="best", bbox_to_anchor=(1.3,1.1), fontsize=6)
        plt.xlabel('time (min)')
        plt.ylabel('Intensity')
        plt.xlim(tmin,tmax)
        plt.ylim(0,int(y[-1:]))
        
        plt.savefig("/Users/bab226/Documents/yale_research/iapp/tht_data/Figures/" + savetitle + "_tht_fit.pdf", dpi=300, transparent=False, bbox_inches='tight')
    
        idx = label
        
        return tlag, t_50, kapp, max_int, idx