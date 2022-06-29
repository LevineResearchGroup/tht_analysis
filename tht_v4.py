#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 10:14:16 2021

@author: Bryan Bogin
Purpose: Plot ThT data from plate reader excell formatted data.

First loads data from excell file. Then, saves data in pandas.
Data is manipulated so time is in minutes.
Data is fit to a sigmoidal curve with 7 paramaters.
The data is normalized by t_50 and by amplitudes so that curve
shape and kinetics can be compared.

"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import cm
import math
import lmfit
import seaborn as sns
import scipy

def get_cmap(n, name='viridis'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

#Body 

#Data Path:
var = pd.read_excel('/Users/bab226/Documents/yale_research/iapp/tht_data/tht_sds_25c_06_26_22.xlsx')

time = var['Time']

#Convert time clock into minutes
time_conv = pd.to_timedelta(time.astype(str), unit='s')
minutes = time_conv.dt.total_seconds()/60

time = minutes

col = var.columns

#Build model
def sigmoidal_curve(time, yi, yf, mi, mf, t_50, tau):
    return(yi+mi*time+((yf+mf*time)/(1+math.e**(-(time-(t_50)/tau)))))

#Fit data
def sigfit(t, y):
    '''
    t = time
    y = signal
    '''
    yi_guess = np.mean(y[:5])
    yf_guess = np.mean(y[len(y)-5:len(y)])
    y_50 = abs(y - yf_guess/2)
    ind = np.argmin(y_50)
    t_50_guess = t[ind]    
    print('Guessing t_50 = %f min' %(t_50_guess))
    
    mi_guess = 0
    mf = 0
    tau_guess = 1
    
    model = lmfit.Model(sigmoidal_curve)
    params = model.make_params()
    params['yi'].set(min=0, value=yi_guess)
    params['yf'].set(min=0, value=yf_guess)
    params['t_50'].set(min=0.01, max=120, value=t_50_guess)
    params['mi'].set(min=0, value=mi_guess)
    params['mf'].set(min=0, value=mf)
    params['tau'].set(min=0, max=5, value=tau_guess)
    
    #weights = 1/err
        
    fitres = model.fit(y, time=t, params=params, method='least_squares')
    # fitres = model.fit(y, t, params=params, method='least_squares',
    #                    weights=weights[t<t_max])
    print('\nList of fitted parameters for %s: \n' % model.name)
    fitres.params.pretty_print(colwidth=10, columns=['value', 'stderr', 'min', 'max'])
    
    print(fitres.fit_report())
    
    return fitres
#%%
#Select Raw Data:

num_reps = 3

rep1 = [col for col in var.columns if 'I' in col]
rep2 = [col for col in var.columns if 'K' in col]
rep3 = [col for col in var.columns if 'M' in col]

headers = []
headers = np.append(headers, rep1[1:7])
headers = np.append(headers, rep2[1:7])
headers = np.append(headers, rep3[1:7])

num_samples = int(len(headers)/num_reps)

num_cols = len(headers)

#alpha = np.linspace(0.7,0.5,len(new_cols))
color = sns.color_palette('deep', num_cols, desat=1)
    
#Title to describe experiment/sample
savetitle = "iapp"

'''
# width = 3.42
# fig, ax = plt.subplots(1, 1, figsize=(width,width/1.2))

# t_50_array = []
    
# Select columns to plot.
for i in range(0, len(headers)):

    minVal = np.amin(var[headers[i]])
    
    y = var[headers[i]] - minVal
    
    #Plotting raw data:
    ax.plot(time, var[headers[i]], color=color[i], label=headers[i])
    ax.set_xlabel('time (min)')
    ax.set_ylabel('Relative Intensity')
    ax.set_xlim(0,1000)
    ax.legend(loc="best", bbox_to_anchor=(1.0,0.8), fontsize=12)
    
    #Sigmoidal fit:
    fitres = sigfit(time, y)

    t_50 = fitres.values["t_50"]
    t_50_array = np.append(t_50_array, fitres.values['t_50'])
    
    yi = fitres.values["yi"]
    yf = fitres.values["yf"]
    '''
    #Plotting transformed data:
    # ax.plot(time, y, color=color[i], label=headers[i])
    # ax.set_xlabel('time (min)')
    # ax.set_ylabel('Relative Intensity')
    # ax.set_xlim(0,1000)
    #ax.legend(loc="best", bbox_to_anchor=(1.0,0.8), fontsize=12)
    
    # ax.plot(time, fitres.best_fit, color = 'm', label = 'fit', linestyle = '--', zorder = 10)

# plt.savefig("/Users/bab226/Documents/yale_research/iapp/tht_data/Figures/" + savetitle + ".png", dpi=300, transparent=False, bbox_inches='tight')


    
#%%
#Plot t50s

# width = 3.42
# fig, ax2 = plt.subplots(1, 1, figsize=(width,width/1.2))

# avg_t_50_array = []
# for i in range(0,num_samples):
#     t50_s1 = t_50_array[i::num_samples]
#     mean = np.mean(t50_s1)
#     avg_t_50_array = np.append(avg_t_50_array, mean)

# ax2.plot(headers[0:num_samples], avg_t_50_array, color = 'm', label=headers[0:num_samples], linestyle = '-')

#%% Average Data

''' 
Set initial ThT intensity to 0
Divide signal by maximum value yf
Then average time points between replicates
Plot averaged data.
'''
signal_norm_avg_array = pd.DataFrame(time)

labels = []
plt_labels = np.array(['35uM', '17uM', '9uM', '4uM', '2uM', '1uM'])
for i in range(0,num_samples):
    samples = headers[i::num_samples]
    
    y = var[samples]
    #Plotting raw data:
    # ax.plot(time, var[headers[i]], color=color[i], label=headers[i])
    # ax.set_xlabel('time (min)')
    # ax.set_ylabel('Relative Intensity')
    # ax.set_xlim(0,1000)
    # ax.legend(loc="best", bbox_to_anchor=(1.0,0.8), fontsize=12)
    
    #Functionalize later?
    signal_norm_array = pd.DataFrame(time)
    for j in range(0, len(samples)):
        
        signal = y[samples[j]]
        
        yi_guess = np.mean(signal[:5])
        yf_guess = np.mean(signal[len(y)-5:len(y)])
        y_50 = abs(signal - yf_guess/2)
        ind = np.argmin(y_50)
        t_50_guess = time[ind]    
        #print('Estimate t_50 = %f min' %(t_50_guess))
        
        #Crude normalization
        signal_norm = (signal-yi_guess)
    
        #Plotting transformed data:
        plt.plot(time, signal_norm, color=color[i], label=plt_labels[i])
        plt.xlabel('time (min)')
        plt.ylabel('Relative Intensity')
        plt.xlim(0,1000)
        plt.legend(loc="best", bbox_to_anchor=(1.3,1.1), fontsize=12)
    
        signal_norm_array = signal_norm_array.join(signal_norm)
        
    plt.savefig("/Users/bab226/Documents/yale_research/iapp/tht_data/Figures/" + savetitle + "_tht.png", dpi=300, transparent=False, bbox_inches='tight')

    name = samples[0][1:]
    signal_norm_avg = signal_norm_array[samples].mean(axis=1)  
    signal_norm_avg = signal_norm_avg.to_frame()
    signal_norm_avg.columns = [name]
    
    signal_norm_avg_array = signal_norm_avg_array.join(signal_norm_avg)
    labels = np.append(labels, name)
    
# Averaged plots:
time = signal_norm_avg_array['Time']

yf_array = []
for z in range(0, num_samples):
    label = signal_norm_avg_array.columns[z+1]
    proc_signal = signal_norm_avg_array[label]
    
    yi_post = np.mean(proc_signal[:5])
    yf_post = np.mean(proc_signal[len(proc_signal)-5:len(proc_signal)])
    y_50 = abs(proc_signal - yf_post/2)
    ind = np.argmin(y_50)
    t_50_post = time[ind]    
    print('Estimate for %s t_50 = %f min' %(plt_labels[z], t_50_post))
    
    yf_array.append(yf_post)
    
width = 3.42
fig, ax3 = plt.subplots(1, 1, figsize=(width,width/1.2))
ax3.plot(time, signal_norm_avg_array[labels], label=plt_labels)
ax3.set_xlabel('time (min)')
ax3.set_ylabel('Relative Intensity')
ax3.set_xlim(0,1000)
ax3.legend(loc="best", bbox_to_anchor=(1.0,0.8), fontsize=12)

plt.savefig("/Users/bab226/Documents/yale_research/iapp/tht_data/Figures/" + savetitle + "_tht_avg.png", dpi=300, transparent=False, bbox_inches='tight')

x = labels[1:len(labels)].astype(float)
y = yf_array[1:len(labels)]

# Plot relative yf vs conc.
width = 3.42
fig, ax4 = plt.subplots(1, 1, figsize=(width,width/1.2))
ax4.scatter(x, y, color = 'm', label=plt_labels)
ax4.set_xlabel('IAPP Concentration')
ax4.set_ylabel('Relative Saturation Intensity')

m, b, r, p_value, std_err = scipy.stats.linregress(x, y)
plt.plot(x, m*x+b)

print("r-squared = %s" %(r))
print("Fit: y = %sx + %s" %(m, b))

plt.savefig("/Users/bab226/Documents/yale_research/iapp/tht_data/Figures/" + savetitle + "_yf.png", dpi=300, transparent=False, bbox_inches='tight')

    
#%% Average Data

'''
NORMALIZED
''' 
'''
Set initial ThT intensity to 0
Divide signal by maximum value yf
Then average time points between replicates
Plot averaged data.
'''
signal_norm_avg_array = pd.DataFrame(time)
plt_labels = np.array(['35uM', '17uM', '9uM', '4uM', '2uM', '1uM'])

labels = []
for i in range(0,num_samples):
    samples = headers[i::num_samples]
    
    y = var[samples]
    #Plotting raw data:
    # ax.plot(time, var[headers[i]], color=color[i], label=headers[i])
    # ax.set_xlabel('time (min)')
    # ax.set_ylabel('Relative Intensity')
    # ax.set_xlim(0,1000)
    # ax.legend(loc="best", bbox_to_anchor=(1.0,0.8), fontsize=12)
    
    #Functionalize later?
    signal_norm_array = pd.DataFrame(time)
    for j in range(0, len(samples)):
        
        signal = y[samples[j]]
        
        yi_guess = np.mean(signal[:5])
        yf_guess = np.mean(signal[len(y)-5:len(y)])
        y_50 = abs(signal - yf_guess/2)
        ind = np.argmin(y_50)
        t_50_guess = time[ind]    
        #print('Estimate t_50 = %f min' %(t_50_guess))
        
        #Crude normalization
        signal_norm = (signal-yi_guess)/(yf_guess-yi_guess)
    
        #Plotting transformed data:
        plt.plot(time, signal_norm, color=color[i], label=plt_labels[i])
        plt.xlabel('time (min)')
        plt.ylabel('Relative Intensity')
        plt.xlim(0,1000)
        plt.legend(loc="best", bbox_to_anchor=(1.3,1.1), fontsize=12)
    
        signal_norm_array = signal_norm_array.join(signal_norm)
        
    plt.savefig("/Users/bab226/Documents/yale_research/iapp/tht_data/Figures/" + savetitle + "_norm_tht.png", dpi=300, transparent=False, bbox_inches='tight')

    name = samples[0][1:]
    signal_norm_avg = signal_norm_array[samples].mean(axis=1)  
    signal_norm_avg = signal_norm_avg.to_frame()
    signal_norm_avg.columns = [name]
    
    signal_norm_avg_array = signal_norm_avg_array.join(signal_norm_avg)
    labels = np.append(labels, name)
    
# Averaged plots:
time = signal_norm_avg_array['Time']

for z in range(0, num_samples):
    label = signal_norm_avg_array.columns[z+1]
    proc_signal = signal_norm_avg_array[label]
    
    yi_post = np.mean(proc_signal[:5])
    yf_post = np.mean(proc_signal[len(proc_signal)-5:len(proc_signal)])
    y_50 = abs(proc_signal - yf_post/2)
    ind = np.argmin(y_50)
    t_50_post = time[ind]    
    print('Estimate for %s t_50 = %f min' %(plt_labels[z], t_50_post))


width = 3.42
fig, ax5 = plt.subplots(1, 1, figsize=(width,width/1.2))
ax5.plot(time, signal_norm_avg_array[labels[0:5]], label=plt_labels[0:5])
ax5.set_xlabel('time (min)')
ax5.set_ylabel('Relative Intensity')
ax5.set_xlim(0,1000)
ax5.legend(loc="best", bbox_to_anchor=(1.0,0.8), fontsize=12)

plt.savefig("/Users/bab226/Documents/yale_research/iapp/tht_data/Figures/" + savetitle + "_norm_tht_avg.png", dpi=300, transparent=False, bbox_inches='tight')

#%%
# n = 0

# t_50_array = []

# #Normalized plot PBS
# plt.figure()
# savetitle = "tht_pbs"
# n = 0
# for cols in new_cols:
    
#     plt.plot(var['Time']/t_50, (var[cols]-yi)/yf, color=color[n])
#     plt.plot(time/t_50, (fitres.best_fit-yi)/yf, color = 'm', label = 'fit', linestyle = '--', zorder = 10)
    
#     #plt.ylim(0.65,0.90)
#     #plt.xlim()
#     plt.xlabel(r't/$t_{50}$')
#     plt.ylabel('Normalized Intensity')
#     #plt.xlim(0,20000)
#     plt.legend(loc="lower right", bbox_to_anchor=(1.0,0.4), fontsize=12)
        
#     n = n + 1
    
# plt.savefig("./Figures/tht_v3/" + savetitle + ".pdf", dpi=300, transparent=False, bbox_inches='tight')


# Make one fig with all plots
# fig, ax = plt.subplots(23, 1, figsize=(width,(width/1.62)*25), sharex=False)

# for i in columns[2:]:
#     plt.subplots_adjust(hspace=0.1)
    
#     ax[n].plot(var['Time'], var[i], color=color(n*50))
    
#     ax[n].set_ylim(0,2500)

#     ax[n].set_title(i)

#     ax[n].set_xlabel('time (min)', labelpad=5)
    
#     ax[n].set_ylabel('A.U.', labelpad=5)
    
#     n = n + 1

#print (varNew)

#varNew.to_csv('newData.csv', index = False)

#plt.plot(varNew[2,:])
#varNew.head()

#os.mkdir("./Figures/")

#plt.savefig("./Figures/" + savetitle + ".pdf", dpi=300, transparent=False, bbox_inches='tight')

