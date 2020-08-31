#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 10:27:46 2019

@author: misiak
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from data_classes import Data_mmr3
from model_classes import Model_ruo2
from ruo2_equation import mmr3_eq_temperature, mmr3_eq_resistance

plt.close('all')
plt.rcParams['text.usetex']=True
def sanitize(string):
    return string.replace('_', ' ')

style = [pe.Normal(), pe.withStroke(foreground='k', linewidth=3)]

data_paths = ('/home/misiak/Data/data_run59/mmr3/MACRT_2019-09-16.csv', )

data = Data_mmr3(data_paths, version='new')

ruo2_labels = ['RuO2 IP2I 1_Meas', 'RuO2 IP2I 2_Meas', 'RuO2 IP2I 3_Meas']

model = Model_ruo2()

param_dict = dict()
for rlab in ruo2_labels:
    param_dict[rlab] = dict()
    
# RuO2 IP2I 1_Meas
param_dict['RuO2 IP2I 1_Meas']['inf40'] = [0.1142, 2843, 2.027]
param_dict['RuO2 IP2I 1_Meas']['sup40'] = [0.4588, 1398, 2.808]
# RuO2 IP2I 2_Meas
param_dict['RuO2 IP2I 2_Meas']['inf40'] = [0.1163, 2744, 2.095]
param_dict['RuO2 IP2I 2_Meas']['sup40'] = [0.5047, 1329, 2.916]
# RuO2 IP2I 3_Meas
param_dict['RuO2 IP2I 3_Meas']['inf40'] = [0.1801, 2181, 2.460]
param_dict['RuO2 IP2I 3_Meas']['sup40'] = [0.6170, 1229, 3.107]


def function_ruo2(ruo2_key, res_array):
    
    # inf to 40mK
    pinf = param_dict[ruo2_key]['inf40']
    temp_inf = mmr3_eq_temperature(pinf, res_array)
    
    # sup to 40mK
    psup = param_dict[ruo2_key]['sup40']
    temp_sup = mmr3_eq_temperature(psup, res_array)
    
    temp_array = np.where(temp_inf<0.040, temp_inf, temp_sup)
    
    return temp_array
 


def function_ruo2_inv(ruo2_key, temp_array):
    
    # inf to 40mK
    pinf = param_dict[ruo2_key]['inf40']
    res_inf = mmr3_eq_resistance(pinf, temp_array)
    
    # sup to 40mK
    psup = param_dict[ruo2_key]['sup40']
    res_sup = mmr3_eq_resistance(psup, temp_array)
    
    res_array = np.where(temp_array<0.04, res_inf, res_sup)
    
    return res_array
 

   
temp_ruo2 = {rlab: function_ruo2(rlab, data.data_dict[rlab]) for rlab in ruo2_labels}

diff_ruo2 = {rlab: np.abs(data.temperature-temp_ruo2[rlab])/data.temperature for rlab in ruo2_labels}

# =============================================================================
# PLOT ruo2 char
# =============================================================================
fig = plt.figure('Plot ruo2 char')
ax = fig.subplots()

for rlab in ruo2_labels:
    
    res_array = np.array(data.data_dict[rlab])
    temp_array = np.array(data.temperature)
    
    temp_cover = np.linspace(temp_array.min(), temp_array.max(), 1000)
    
    line, = ax.plot(
            temp_array, res_array,
            label=sanitize(rlab),
            ls='none', marker='+', markersize=12, alpha=0.3
    )
    
    ax.plot(
            temp_cover, function_ruo2_inv(rlab, temp_cover),
            label=sanitize(rlab),
            ls='-', color=line.get_color(),
            path_effects=style,
            zorder=10,    
    )
    
ax.grid(True)
ax.set_yscale('log')
ax.set_xscale('log')
ax.legend()

ax.set_title('RuO2 Law Calibration R=f(T)')
ax.set_xlabel('Temperature [K]')
ax.set_ylabel('Resistance [$\Omega$]')
fig.tight_layout()





# =============================================================================
# PLOT vs TIME
# =============================================================================
fig = plt.figure('Plot vs Time')
axes = fig.subplots(nrows=2, sharex=True)

axes[0].plot(data.time, data.temperature, label='MC')

axes[1].axhline(1, ls='-', color='k')
axes[1].axhline(0.01, ls='-.', color='k')

for rlab in ruo2_labels:
    line, = axes[0].plot(data.time, temp_ruo2[rlab], label=sanitize(rlab))

    axes[1].plot(data.time, diff_ruo2[rlab], label=sanitize(rlab), color=line.get_color())

for ax in axes:
    ax.set_yscale('log')
    ax.grid(True)
    ax.legend()

axes[0].set_title('RuO2 Calibration vs Time')
axes[0].set_ylabel('Temperature [K]')
axes[1].set_ylabel('Relative Residual dT/T')
axes[-1].set_xlabel('Time Unix [s?]')

fig.tight_layout()
fig.subplots_adjust(hspace=0)

# =============================================================================
# PLOT vs Temperature
# =============================================================================
fig = plt.figure('Plot vs Temperature')
axes = fig.subplots(nrows=2, sharex=True)

axes[0].plot(
        data.temperature, data.temperature,
        label='MC', ls='none', marker='o',
)

axes[1].axhline(1, ls='-', color='k')
axes[1].axhline(0.01, ls='-.', color='k')

for rlab in ruo2_labels:
    line, = axes[0].plot(
            data.temperature, temp_ruo2[rlab],
            label=sanitize(rlab),
            ls='none', marker='+',
    )
    
    axes[1].plot(
            data.temperature, diff_ruo2[rlab],
            label=sanitize(rlab),
            ls='none', marker='+', color=line.get_color()
    )

for ax in axes:
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True)
    ax.legend()

axes[0].set_title('RuO2 Calibration vs Temperature')
axes[0].set_ylabel('Temperature [K]')
axes[1].set_ylabel('Relative Residual dT/T')
axes[-1].set_xlabel('Temperature MC [K]')

fig.tight_layout()
fig.subplots_adjust(hspace=0)
