#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 11:23:49 2022

@author: bab226

Purpose:
Anaylsis of ThT curves.
Options for traditional sigmoid fitting and alternative ad hoc fittings.

Determins kinetic paramaters such as t50, tlag, kapp, and max intensity:
    
t50 = time to half-max signal
tlag = time delay before exponential rise defined by intercept of lag and exp
phase.
kapp = observed kinetic rate of aggregation
max intensity = highest observed signal

"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import cm
import math
import lmfit
import seaborn as sns
from scipy.stats import linregress
from uncertainties import ufloat
import tht_fitting as fitme
import numdifftools as nd


"""""""""""""""
""""""BODY""""""""
"""""""""""""""

'''Save Directory'''
savedir = "/Users/bab226/Documents/yale_research/iapp/tht_data/"

savetitle = "iapp_adm16_07_03_22"

'''Number of replicates'''
num_reps = 3

'''Time of assay in minutes'''
mes_time = 360

"""Import excell file as pandas df"""
df = pd.read_excel('/Users/bab226/Documents/yale_research/iapp/tht_data/tht_iapp_30uM_screen_06_29_22.xlsx')

headers = df.columns

'''Define Time Period'''
tmin = 0
tmax = 360

'''Color MAP'''
color = sns.color_palette('deep', len(headers), desat=1)

'''Select corresponding index range'''
#Convert time clock into minutes
time_conv = pd.to_timedelta(df['Time'].astype(str), unit='s')
time = time_conv.dt.total_seconds()/60

interval = int(np.round(time[1]-time[0]))

ind_min = int(tmin/interval)
ind_max = int(tmax/interval)

'''Define number of samples'''
num_samples = int(len(df.columns[2:])/num_reps)

#key = [col for col in df.columns if 'pbs' in col]
plt.figure()
'''Raw Data Plot'''
z = 0
step = 0
plt.figure()
for col in df.columns[2:]:

    color_ind = list(range(0,len(color), 3))
    x = time[ind_min:ind_max]
    y = df[col][ind_min:ind_max]
    label = col
    step = int(z/3)
    colors = color[color_ind[step]]
        
    plt.plot(x, y, color=colors, label=label)
    plt.legend(loc="best", bbox_to_anchor=(1.0,1.1), fontsize=12)
    plt.xlabel('time (min)', fontsize=16)
    plt.ylabel('Intensity', fontsize=16)
    plt.xlim(tmin,tmax)
    plt.ylim(0,10000)

    z = z + 1

plt.savefig("/Users/bab226/Documents/yale_research/iapp/tht_data/Figures/" + savetitle + "_tht.pdf", dpi=300, transparent=False, bbox_inches='tight')
    
'''Fit ThT Curves'''
tlag_arr = []
t50_arr = []
tau_arr = []
id_arr = []
vt50_arr = []
vt50_norm_arr = []

good_set = []
# tlag_arr_err = []
# t50_arr_err = []
# kapp_arr_err = []
# max_int_arr_err = []
j = 0
'''Fitting and analysis'''
plt.figure()
#key = [col for col in df.columns if 'pbs' in col]
for col in df.columns[2:]:
    x = time[ind_min:ind_max]
    y = df[col][ind_min:ind_max]
    label = col
    
    t50_guess, tau_guess, vt50_guess, c = fitme.guess_fit(x, y, label, interval, savetitle, headers, j)
    
    ###***FIXME***### Fitting with sigmoidal curves were poor
    #fitres = fitme.sigfit(x, y, label, savetitle, i, tau_guess, c)
        
    if str(t50_guess) != 'nan':
        
        t50 = t50_guess
        tau = tau_guess
        vt50 = vt50_guess
        tlag = t50 - 2*tau
        
    #     t50_guess, tau_guess, c = fitme.guess_fit(x, y, label, savetitle, headers, i, t50)
       
        t50_arr = np.append(t50_arr, t50)
        tau_arr = np.append(tau_arr, tau)
        tlag_arr = np.append(tlag_arr, tlag)
        vt50_arr = np.append(vt50_arr, vt50)
        id_arr = np.append(id_arr, label)
    
        # '''Normalization'''
        # if y[int(t50)] > 0:
        #     good_set = np.append(good_set, label)
            
        #     t50 = int(t50)
        #     vt50_norm = (vt50/y[int(t50)])*t50
            
        # else: 
        #     vt50_norm = np.nan
            
        # vt50_norm_arr = np.append(vt50_norm_arr, vt50_norm)
        
        j = j + 1
        
    else:
        print ("Is NAN. Skip ME!")
        
plt.savefig("/Users/bab226/Documents/yale_research/iapp/tht_data/Figures/" + savetitle + "_tht_fit.pdf", dpi=300, transparent=False, bbox_inches='tight')
            
plt.figure()
'''Plot lag half-normalized graph'''
for i in range(0, len(vt50_norm_arr)):
    if y[int(t50_arr[i])] > 0:
        col = good_set[i]
        
        x = time[ind_min:ind_max]
        y = df[col][ind_min:ind_max]
    
        xnorm = x - tlag_arr[i]
        ynorm = y/y[int(t50_arr[i])]
            
        t50_int = fitme.find_nearest(x, t50)
        
        #Color map
        color = sns.color_palette('deep', len(vt50_norm_arr), desat=1)
        color_ind = list(range(0,len(color), 3))
        step = int(i/3)
        color = color[step]
        
        plt.plot(xnorm, ynorm, label=col, color=color)
        plt.legend(loc="best", bbox_to_anchor=(1.0,1.0), fontsize=12)
        plt.xlabel('t-tlag (min)', fontsize=16)
        plt.ylabel('y/y50', fontsize=16)
        plt.xlim(0, tmax - tlag_arr[i])
        plt.ylim(0,2.2)


plt.savefig("/Users/bab226/Documents/yale_research/iapp/tht_data/Figures/" + savetitle + "_tht_half_norm.pdf", dpi=300, transparent=False, bbox_inches='tight')

plt.figure()
'''Plot normalized graph'''
for i in range(0, len(vt50_norm_arr)):
    if y[int(t50_arr[i])] > 0:
        col = good_set[i]
        
        x = time[ind_min:ind_max]
        y = df[col][ind_min:ind_max]
    
        xnorm = x/t50_arr[i]
        ynorm = y/y[int(t50_arr[i])]
            
        t50_int = fitme.find_nearest(x, t50)
        
        #Color map
        color = sns.color_palette('deep', len(vt50_norm_arr), desat=1)
        color_ind = list(range(0,len(color), 3))
        step = int(i/3)
        color = color[step]
        
        plt.plot(xnorm, ynorm, label=col, color=color)
        plt.legend(loc="best", bbox_to_anchor=(1.0,1.0), fontsize=12)
        plt.xlabel('t/t50 (min)', fontsize=16)
        plt.ylabel('y/y50', fontsize=16)
        plt.xlim(tmin/x[t50],2.2)
        plt.ylim(0,2.0)

plt.savefig("/Users/bab226/Documents/yale_research/iapp/tht_data/Figures/" + savetitle + "_tht_norm.pdf", dpi=300, transparent=False, bbox_inches='tight')
    

'''Summary'''
# tlag_arr = [x for x in tlag_arr if str(x) != 'nan']
# t_50_arr = [x for x in t_50_arr if str(x) != 'nan']
# kapp_arr = [x for x in kapp_arr if str(x) != 'nan']
# max_int_arr = [x for x in max_int_arr if str(x) != 'nan']
# id_arr = [x for x in id_arr if str(x) != 'nan']
'''Averaged data'''

'''Define lists'''


tlag_all = []
t50_all = []
kapp_all = []
max_int_all = []
vt50_all = []

tlag_all_err = []
t50_all_err = []
kapp_all_err = []
max_int_all_err = []
vt50_all_err = []


id_base_arr = []
'''Save all averaged fitting data to list'''


for i in range(0, num_reps*int(len(id_arr)/num_reps), num_reps):
    
    id_base = id_arr[i][:-1]
    id_base_arr = np.append(id_base_arr, id_base)
    
    tlag_avg = np.mean(tlag_arr[i:i+3])
    t50_avg = np.mean(t50_arr[i:i+3])
    vt50_avg = np.mean(vt50_arr[i:i+3])
    
    tlag_all = np.append(tlag_all, tlag_avg)
    t50_all = np.append(t50_all, t50_avg)
    vt50_all = np.append(vt50_all, vt50_avg)
    
    tlag_err = np.std(tlag_arr[i:i+3])
    t50_err = np.std(t50_arr[i:i+3])
    vt50_err = np.std(vt50_arr[i:i+3])

    tlag_all_err = np.append(tlag_all_err, tlag_err)
    t50_all_err = np.append(t50_all_err, t50_err)
    vt50_all_err = np.append(vt50_all_err, vt50_err)
    
'''Write lists to file'''
# with open(savedir + 'tlag' + '.txt', "w" ) as f:
#     f.write('ID \t tlag \t err \n')
#     for ID,R,err in zip(id_arr, tlag_arr, tlag_arr_err):
#         f.write('%s \t %.3f \t %.3f \n' %(ID,R,err))
   
# with open(savedir + 't50' + '.txt', "w" ) as f:
#     f.write('ID \t t50 \t err \n')
#     for ID,R,err in zip(id_arr, t50_arr, t50_arr_err):
#         f.write('%s \t %.3f \t %.3f \n' %(ID,R,err))

with open(savedir + 'tlag_avg' + '.txt', "w" ) as f:
    f.write('ID \t tlag \t err \n')
    for ID,R,err in zip(id_base_arr, tlag_all, tlag_all_err):
        f.write('%s \t %.3f \t %.3f \n' %(ID,R,err))
        
with open(savedir + 't50_avg' + '.txt', "w" ) as f:
    f.write('ID \t t50 \t err \n')
    for ID,R,err in zip(id_base_arr, t50_all, t50_all_err):
        f.write('%s \t %.3f \t %.3f \n' %(ID,R,err))

with open(savedir + 'vt50_avg' + '.txt', "w" ) as f:
    f.write('ID \t vt50 \t err \n')
    for ID,R,err in zip(id_base_arr, vt50_all, vt50_all_err):
        f.write('%s \t %.3f \t %.3f \n' %(ID,R,err))

'''
plt.figure()
#Plot lag half-normalized graph avg
yavg_arr = pd.DataFrame(time)
for i in range(0, num_reps*(len(id_base_arr)), num_reps):
    col = df.columns[2+i:2+i+num_reps]   #skip time and temp. columns
    x = time[ind_min:ind_max]
    y = df[col][ind_min:ind_max]
    label = id_base_arr[int(i/3)]
   
    yavg = y.mean(axis=1)
    yavg = yavg.to_frame()
    yavg.columns = [label]
    
    yavg_arr = yavg_arr.join(yavg)
    
    #plot 
    xnorm = x - tlag_all[int(i/3)]
    ynorm = yavg[label]/yavg[label][int(t50_all[int(i/3)])]
    
    #Color map
    color = sns.color_palette('deep', len(vt50_norm_arr), desat=1)
    color_ind = list(range(0,len(color), 3))
    step = int(i/3)
    color = color[step]
    
    plt.plot(xnorm, ynorm, label=label, color=color)
    plt.legend(loc="best", bbox_to_anchor=(1.0,1.0), fontsize=12)
    plt.xlabel('t-tlag (min)', fontsize=16)
    plt.ylabel('y/y50', fontsize=16)
    plt.xlim(0, tmax)
    plt.ylim(0,2.2)

plt.savefig("/Users/bab226/Documents/yale_research/iapp/tht_data/Figures/" + savetitle + "_tht_half_norm_avg.pdf", dpi=300, transparent=False, bbox_inches='tight')
    
plt.figure()
#Plot normalized graph avg
yavg_arr = pd.DataFrame(time)
for i in range(0, num_reps*(len(id_base_arr)), num_reps):
    col = df.columns[2+i:2+i+num_reps]   #skip time and temp. columns
    x = time[ind_min:ind_max]
    y = df[col][ind_min:ind_max]
    label = id_base_arr[int(i/3)]
   
    yavg = y.mean(axis=1)
    yavg = yavg.to_frame()
    yavg.columns = [label]
    
    yavg_arr = yavg_arr.join(yavg)
    
    #plot 
    xnorm = x/t50_all[int(i/3)]
    ynorm = yavg[label]/yavg[label][int(t50_all[int(i/3)])]
    
    #Color map
    color = sns.color_palette('deep', len(vt50_norm_arr), desat=1)
    color_ind = list(range(0,len(color), 3))
    step = int(i/3)
    color = color[step]
    
    plt.plot(xnorm, ynorm, label=label, color=color)
    plt.legend(loc="best", bbox_to_anchor=(1.0,1.0), fontsize=12)
    plt.xlabel('t-tlag (min)', fontsize=16)
    plt.ylabel('y/y50', fontsize=16)
    plt.xlim(0, 2.2)
    plt.ylim(0,2.2)

plt.savefig("/Users/bab226/Documents/yale_research/iapp/tht_data/Figures/" + savetitle + "_tht_norm_avg.pdf", dpi=300, transparent=False, bbox_inches='tight')
    
'''











