#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: misiak

"""

import matplotlib.pyplot as plt

import scipy.optimize as op
import numpy as np
import mcmc_red as mcr

from data_classes import Data_mmr3
from model_classes import Model_ruo2
from fitting_class import Manual_fitting

plt.close('all')
plt.rcParams['text.usetex']=True
def sanitize(string):
    return string.replace('_', ' ')


#    data_paths = ('test/MACRT_2019-07-15.csv', 'test/MACRT_2019-07-16.csv')
data_paths = ('/home/misiak/Data/data_run59/mmr3/MACRT_2019-09-16.csv', )

data_ruo2 = Data_mmr3(data_paths, version='new')

ruo2_labels = ['RuO2 IP2I 1_Meas', 'RuO2 IP2I 2_Meas', 'RuO2 IP2I 3_Meas']

#%%
# =============================================================================
# TIME CUT
# =============================================================================
t_inf = 3.65146e9 +3500
t_sup = 3.6514e9 +90000

time_cut = np.logical_and( t_inf < data_ruo2.time, data_ruo2.time < t_sup)


# TIME CUT
fig, axes = plt.subplots(nrows=2, sharex=True, num='time cut plot')

axes[0].plot(
        data_ruo2.time,
        data_ruo2.temperature,
        ls='none', marker='+', color='r'
)

axes[0].plot(
        data_ruo2.time[time_cut],
        data_ruo2.temperature[time_cut],
        ls='none', marker='+', color='r'
)

for rlab in ruo2_labels:
    res_array = data_ruo2.data_dict[rlab]
    axes[1].plot(
            data_ruo2.time,
            res_array,
            label=sanitize(rlab)
    )

for ax in axes:
    ax.set_yscale('log')
    ax.axvline(t_inf, color='k')
    ax.axvline(t_sup, color='k')
    ax.legend()
    ax.grid()

#%%
# =============================================================================
# ROI CUT
# =============================================================================
fig, ax = plt.subplots(nrows=3, num='roi cut plot', sharex=True)

temp_array = data_ruo2.temperature[time_cut]
x_indexes = np.arange(len(temp_array))

ax[0].plot(
        x_indexes,
        temp_array,
        ls='none', marker='+', color='r', 
        label='all'
)

grad_temp = np.gradient(data_ruo2.temperature[time_cut])
pseudo_grad =  np.abs(grad_temp) / data_ruo2.temperature[time_cut]

ax[1].plot(
        x_indexes,
        pseudo_grad,
        ls='none', marker='+', color='b',
        label='all'
)

grad_cut = pseudo_grad < 5e-3

grad_indexes = x_indexes[grad_cut] - np.roll(x_indexes[grad_cut], shift=1)
x_starts = np.where(grad_indexes != 1)[0]

section_indexes = [sec for sec in np.split(x_indexes[grad_cut], x_starts) if len(sec)!=0]
section_temp = [temp_array[sec] for sec in section_indexes]

for sec, tem in zip(section_indexes, section_temp):
    line_ax0 = ax[0].plot(
            sec,
            tem,
            ls='-', color='b'
    )

    line_ax1 = ax[1].plot(
            sec,
            pseudo_grad[sec],
            ls='-', color='r'
    )

line_ax0[0].set_label('sec cut')
line_ax1[0].set_label('sec cut')

buffer = 5
roi_indexes = [sec[buffer:] for sec in section_indexes if len(sec)>5]
roi_temp = [sec[buffer:] for sec in section_temp if len(sec)>5]

for sec, tem in zip(roi_indexes, roi_temp):
    line_ax0 = ax[0].plot(
            sec,
            tem,
            ls='-', color='lightskyblue'
    )
line_ax0[0].set_label('roi cut')

roi_res_dict = dict()
for rlab in ruo2_labels:
    res_array = data_ruo2.data_dict[rlab][time_cut]
    
    line = ax[2].plot(x_indexes, res_array,
               label=sanitize(rlab) + ' time cut', alpha=0.1)[0]
    
    color = line.get_color()
    
    for sec in section_indexes:
        line = ax[2].plot(
                sec, res_array[sec],
                alpha=0.3, color=color
        )[0]
    line.set_label(sanitize(rlab) + ' section cut')
    
    for roi in roi_indexes:
        line = ax[2].plot(
                roi, res_array[roi],
                marker='+', color=color
        )[0]
    line.set_label(sanitize(rlab) + ' roi')    

    roi_res_dict[rlab] = [res_array[roi] for roi in roi_indexes]

for a in ax:
    a.grid(True)
    a.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
    a.set_yscale('log')
fig.tight_layout()


#%%

# mean and std on the steps
mean_temp = [np.mean(tem) for tem in roi_temp]
mean_res_dict = dict()
std_res_dict = dict()
len_res_dict = dict()
for rlab in ruo2_labels:
    mean_res_dict[rlab] = [np.mean(res) for res in roi_res_dict[rlab]]
    std_res_dict[rlab] = [np.std(res)/len(res)**0.5 for res in roi_res_dict[rlab]]

# =============================================================================
# EXPO PLOT
# =============================================================================
rlab = ruo2_labels[0]
data_array = mean_res_dict[rlab]
std_array = std_res_dict[rlab]

# used for the fitting
cut = np.array(mean_temp) < 45e-3
#cut = np.array(mean_temp) > 35e-3
mean_temp = np.array(mean_temp)[cut]
data_array = np.array(data_array)[cut]
std_array = np.array(std_array)[cut]

cover_temp = np.linspace(np.min(mean_temp), np.max(mean_temp), 1000)

model = Model_ruo2()

def function_inv(param, res_array):
    a1, a2, a3 = param
    return a1 / ( np.log( res_array/a2 ) )**a3

def model_fun(param):
    term0 = model.function(param, cover_temp)
    delta_t = mean_temp - function_inv(param, data_array)
    term1 = delta_t / mean_temp
#    term1 = (data_array - model.function(param, mean_temp))/data_array
    return [term0, term1]
   
model_fun_array = model_fun(model.parameters_0)

def chi2_ruo2(param):
    mod_array = model.function(param, mean_temp)
    chi2_ruo2 = mcr.chi2_simple(mod_array, data_array, std_array)
    return chi2_ruo2


fig = plt.figure('ruo2 fitting plot')
axes = fig.subplots(nrows=2, sharex=True)

axes[0].plot(mean_temp, data_array, ls='none', marker='+', label='data')
line0, = axes[0].plot(
        cover_temp,
        model_fun_array[0],
        label='model'
)
axes[0].set_yscale('log')
axes[0].set_xscale('log')

line1, = axes[1].plot(
        mean_temp,
        model_fun_array[1],
        label='dT/T',
        ls='none', marker='+'
)    

for ax in axes:
    ax.grid(True)
    ax.legend()

def chi2_legend(param):
    chi2 = chi2_ruo2(param)
    dof = data_array.shape[0]
    for l in [line0, line0]:
        l.set_label('model\n$\chi^2_n=${:.3e}'.format(chi2/dof))
        l.axes.legend()

mfit = Manual_fitting([line0,line1], model_fun, model.parameters_0,
                      callback=chi2_legend)

result = op.minimize(chi2_ruo2, model.parameters_0,
                     method='Nelder-Mead',
                     options={'maxiter':1000})

mfit.set_param(result.x)

## =============================================================================
## FINAL RESULTS
## =============================================================================
#model = Model_ruo2()
#
#mean_temp = np.array([np.mean(tem) for tem in roi_temp])
#cut_label = ('12mK to 40mK', '40mK to 300mK')
#cut_truth = (mean_temp<45e-3, mean_temp>35e-3)
#
#for rlab in ruo2_labels:
#    print(rlab)
#    for cut,clab in zip(cut_truth, cut_label):
#        print(clab)
#        temp_array = np.array(mean_temp)[cut]
#        data_array = np.array(mean_res_dict[rlab])[cut]
#        std_array = np.array(std_res_dict[rlab])[cut]
#
#        def chi2_ruo2(param):
#            mod_array = model.function(param, temp_array)
#            chi2_ruo2 = mcr.chi2_simple(mod_array, data_array, std_array)
#            return chi2_ruo2
#    
#        result = op.minimize(chi2_ruo2, model.parameters_0,
#                             method='Nelder-Mead',
#                             options={'maxiter':1000})    
#    
#        assert result.success
#        print(result.x)
#    
#    
#    
    
