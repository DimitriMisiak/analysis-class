#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: misiak

"""
import os
import importlib.util

import numpy as np
import matplotlib.pyplot as plt

import ethem as eth

from data_classes import Data_vi, plot_vi

class Model_vi(object):
    
    data_type = 'VI characteristic'
    model_type = 'Electro-Thermal Modelization'
    
    _flag_ready = False
    conditionnal_parameters = tuple()
    parameters = tuple()
    
    def __init__(self, config_path, data, label):
        
        self.config_path  = config_path
        self.label = ' '.join((self.model_type, label))
        
        # associating data
        self.data = data
        self.temp_list = self.data.temp_list
        
        
        # import the configuration file as the module "config"
        basename = os.path.basename(config_path)
        spec = importlib.util.spec_from_file_location(basename, config_path)
        config = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(config)
        
        self.eval_dict = config.get_eval_dict()
        self.system = config.get_system()
    
    @property
    def all_parameters(self):
        return self.conditionnal_parameters + self.parameters
    
    def lock_model(self):
        self._flag_ready = True
        
        # checking if there is two conditionnal parameters: temperature and current
        assert len(self.conditionnal_parameters) == 2
        
        self.steady_state_fun = eth.solve_sse_param(
                model.system, model.all_parameters, model.eval_dict
        )
    
    
    def function(self, param):
        """ Return the model array from the given set of parameters.
        """
        if self._flag_ready == False:
            raise Exception('Model is not ready. Lock model to proceed.')

        model_iv_list = list()
        
        for i,temp in enumerate(self.temp_list):
    
            i_array, _ = self.data.vi_arrays[i]
            
            model_iv = list()
            for curr in i_array:
                
                p_top = (temp, curr)
                p_full = p_top + tuple(param)
                
                model_iv.append(self.steady_state_fun(p_full).x[-1])
                
            model_iv_array = np.vstack((i_array, model_iv))
            model_iv_list.append(model_iv_array)
    
        return model_iv_list
    
    
    def fake_data(self, param, sigma_coeff=0.1):
        """ Add a gaussian noise depending on the model_data.
        """
        
        model_data = self.function(param)
        
        fake_data = list()
        for md in model_data:
    
            x_data, y_data = md
    
            sigma = y_data * sigma_coeff
    
            blob = np.random.normal(0,sigma,sigma.shape)
    
            noisy_data = y_data + blob
    
            fake_data.append(np.vstack((x_data, noisy_data)))
    
        return fake_data
  
    
if __name__ == '__main__':
    
    plt.close('all')
    plt.rcParams['text.usetex']=True

    data_path = 'test/vi_run57_red80_late.csv'
    data_iv = Data_vi(data_path, 'RED80')
    
    config_path = 'test/config_ethem.py'
    
    model = Model_vi(config_path, data=data_iv, label='RED80')

    syst = model.system
    
    cryo = syst.Thermostat_b
    capa = syst.Capacitor_f
    elntd = syst.Resistor_ntd
    leak = syst.ThermalLink_leak
    epcoup = syst.ThermalLink_ep
    glue = syst.ThermalLink_glue

    # "primary" parameters
    model.conditionnal_parameters = (cryo.temperature, capa.current,)

#    # "secondary" parameters
#    model.parameters = (
#            elntd.R0,
#            elntd.T0,
#            leak.cond_alpha,
#            epcoup.cond_alpha,
#    )

    model.lock_model()
    
    modvi_arrays = model.function([])
    
    fig = plot_vi(model.temp_list, modvi_arrays, None, num=model.data.label)
    
    fake_arrays = model.fake_data([], 0.1)
    
    plot_vi(model.temp_list, fake_arrays, 0.2, num=model.data.label,
            ls='none', marker='.')
    
#    model.data.plot(num=model.data.label, ls='none', marker='.')
    
    plt.grid()
    plt.xlabel('Voltage [V]')
    plt.ylabel('Current [A]')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=2)
    plt.title(fig.get_label())
    plt.tight_layout(rect=(0., 0., 1., 1.))    