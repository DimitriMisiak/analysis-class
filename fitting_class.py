#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: misiak

"""
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox, Button
import numpy as np

import mcmc_red as mcr

from data_classes import Data_mmr3
from model_classes import Model_ruo2


class Manual_fitting():
    
    text_boxes = list()
    lines = list()
    func = lambda x: list()
    paramnow = list()
    paraminit = list()
    paramprevious = list()
    fig = None
    callback = lambda x: None
    lines_previous = list()
    
    def _update(self, val):
        # getting the values in the text boxes
        param = [float(tbox.text) for tbox in self.text_boxes]
        
        # no update if the same parameters as before
        if param == self.paramnow:
            print('No update. Same param={}'.format(param))
            return None
        
        # saving the previous set of parameters
        self.paramprevious = self.paramnow.copy()
        self.paramnow = param
        
        # evaluating the func for the entered parameters
        new_data = self.func(self.paramnow)
        
        # replacing the ydata of the lines
        for l, nd in zip(self.lines, new_data):
            l.set_ydata(nd)
        
        # autoscale of the axes
        for l in self.lines:
            l.axes.relim()
            l.axes.autoscale(True)
            
        # refresh the figure    
        if self.fig:
            self.fig.canvas.draw_idle()
            
        # callback function
        self.callback(param)
        
        # explicit sanity check
        print('Updated with param={}'.format(param))
        
    def _reset(self, event):
        for p0, text_box in zip(self.paraminit, self.text_boxes):
            text_box.set_val(str(p0))
    
    def _previous(self, event):
        # uses an auxiliary list to swap list contents
        param_aux = self.paramprevious.copy()
        self.paramprevious = self.paramnow.copy()
        for ppre, tbox in zip(param_aux, self.text_boxes):
            tbox.set_val(str(ppre))
        # to finish, update, some memoization might be possible here
        self._update(0)
    
    def __init__(self, lines, func, paraminit, callback=None):
        
        self.lines = lines
        self.func = func
        self.fig = lines[0].get_figure()
        self.paraminit = list(paraminit)
        self.nparams = len(paraminit)
        if callable(callback):
            self.callback = callback
        
        # reserving some place for the widgets
        rightl = 0.7
        topl = 0.75
        self.fig.subplots_adjust(right=rightl)
        
        resetax = plt.axes([rightl+0.1, topl, 0.1, 0.045])
        self.reset_button = Button(resetax, 'Reset', hovercolor='0.975')
        self.reset_button.on_clicked(self._reset)   
        
        updatax = plt.axes([rightl+0.1, topl-0.05, 0.1, 0.045])
        self.update_button = Button(updatax, 'Update', hovercolor='0.975')
        self.update_button.on_clicked(self._update)        
        
        for i,p in enumerate(self.paraminit):
            bot_level = topl-0.1-0.05*i
            axbox = self.fig.add_axes([rightl+0.1, bot_level, 0.15, 0.045])
            text_box = TextBox(axbox, 'param{}'.format(i), initial=str(p))    
            (self.text_boxes).append(text_box)

        previax = plt.axes([rightl+0.1, bot_level-0.05, 0.1, 0.045])
        self.previous_button = Button(previax, 'Previous', hovercolor='0.975')
        self.previous_button.on_clicked(self._previous)

        # finish with initial update
        self._update(0)

    def set_param(self, param):
        # set parameters from command prompt
        for p, tbox in zip(param, self.text_boxes):
            tbox.set_val(str(p)) 
        self._update(0)
        
        
if __name__ == '__main__':
    
    plt.close('all')
    plt.rcParams['text.usetex']=True

    data_paths = ('test/MACRT_2019-07-15.csv', 'test/MACRT_2019-07-16.csv')
#    data_paths = ('/home/misiak/Data/data_run59/mmr3/MACRT_2019-09-16.csv', )

    data = Data_mmr3(data_paths, version='old')
    model = Model_ruo2()
    
    data_array = data.data_dict['MMR3-205_1_Meas']
    std_array = data_array * 0.01

    cover_temp = np.linspace(
            np.min(data.temperature),
            np.max(data.temperature),
            1000
            )

    def model_fun(param):
        term0 = model.function(param, cover_temp)
        term1 = (data_array - model.function(param, data.temperature))/data_array
        return [term0, term1]
   
    model_fun_array = model_fun(model.parameters_0)

    def chi2_ruo2(param):
        mod_array = model.function(param, data.temperature)
        chi2_ruo2 = mcr.chi2_simple(mod_array, data_array, std_array)
        return chi2_ruo2
        
    fig = plt.figure()
    axes = fig.subplots(nrows=2, sharex=True)
    
    axes[0].plot(data.temperature, data_array,
                 ls='none', marker='+', label='data')
    
    line0, = axes[0].plot(
            cover_temp,
            model_fun_array[0],
            label='model'
    )
    axes[0].set_yscale('log')
    axes[0].set_xscale('log')
    
    line1, = axes[1].plot(
            data.temperature,
            model_fun_array[1],
            label='relative residual',
            ls='none', marker='+'
    )    
    
    for ax in axes:
        ax.grid(True)
        ax.legend()

    def chi2_legend(param):
        
        chi2 = chi2_ruo2(param)
        dof = data.nsamples
        for ax in axes:
            ax.legend(title='$\chi^2=${:.3e}\ndof={:.3e}'.format(chi2, dof)) 

    mfit = Manual_fitting(
            [line0,line1],
            model_fun,
            model.parameters_0,
            callback=chi2_legend
    )
