#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: misiak

"""

import matplotlib.pyplot as plt

import scipy.optimize as op

import mcmc_red as mcr


from data_classes import Data_vi, plot_vi
from model_classes import Model_vi

class Analysis_vi(object):
    
    def __init__(self, data_path, model_path, label):
        
        self.data = Data_vi(data_path, label, error=0.1)
        self.model = Model_vi(model_path, self.data, label)


    def lock_analysis(self, compared_arrays=None, compared_errors=None):
        
        if not self.model._flag_ready:
            self.model.lock_model()
        
        if compared_arrays is None:
            compared_arrays = self.data.vi_arrays
            
        if compared_errors is None:
            compared_errors = self.data.std_arrays
            
        self.compared_arrays = compared_arrays
        self.compared_errors = compared_errors        
        
        def chi2_fun(param, check_print=True):
            
            chi2 = self.model.comparator(
                    param,
                    data_arrays=self.compared_arrays,
                    data_errors=self.compared_errors,
                    check_print=check_print
            )
            
            return chi2
        
        self.chi2_fun = chi2_fun
   
    def manual_fitting(self, num='Manual Fitting'):
        
        data_arrays = self.compared_arrays
        error = self.compared_errors
            
        param0 = self.model.parameters_0
        
        model_arrays = self.model.function(param0)
        
        temp_list = self.model.temp_list
        
        plot_vi(temp_list, data_arrays, error, num=num,
                mode='errorbar', marker='.', ls='none')
        
        fig = plot_vi(temp_list, model_arrays, error, num=num,
                      mode='plot')
        
        chi2 = ana.chi2_fun(param0)
        
        ax = fig.get_axes()[0]
        
        ax.grid()
        ax.set_xlabel('Voltage [V]')
        ax.set_ylabel('Current [A]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(title='$\chi^2=${:.3e} ; dof={}'.format(chi2, ana.model.data.nsamples),
                   loc='upper left', bbox_to_anchor=(1, 1), ncol=2)
        ax.set_title(fig.get_label ())
        fig.tight_layout(rect=(0., 0., 1., 1.))   
    
        return fig


    def minimizer_fitting(self, param_init=None, check_print=False):
        
        if param_init is None:
            param_init = self.model.parameters_0

        param0 = param_init
        
        result = op.minimize(ana.chi2_fun, param0, args=(check_print,), 
                             method='nelder-mead')
        
        return result
        
    
    def mcmc_fitting(self, nsteps=10, param_init=None,
                     sampler_dir='mcmc_sampler/autosave'):
        
        if param_init is None:
            param_init = self.model.parameters_0
    
        # running the mcmc analysis
        bounds = [(p/1.1, p*1.1) for p in param_init]
        
        sampler = mcr.mcmc_sampler(
                self.chi2_fun,
                bounds, 
                nsteps=nsteps,
                path=sampler_dir
        )
        
        return sampler

if __name__ == '__main__':
    
    plt.close('all')
    plt.rcParams['text.usetex']=True

    data_path = 'test/vi_run57_red80_late.csv'
    
    model_path = 'test/config_ethem.py'
    
    label = 'RED80'
    
    ana = Analysis_vi(data_path, model_path, label)
    
    system = ana.model.system
    
    ana.model.conditionnal_parameters = (
            system.Thermostat_b.temperature,
            system.Capacitor_f.current
    )
    
    ana.model.parameters = (
            system.Resistor_ntd.R0,
            system.Resistor_ntd.T0,
            system.ThermalLink_leak.cond_alpha,
            system.ThermalLink_ep.cond_alpha,
    )

    ana.model.lock_model()
    fake_arrays = ana.model.fake_data(ana.model.parameters_0)
    fake_errors = [a*0.1 for a in fake_arrays]
    
    
    ana.lock_analysis(fake_arrays, fake_errors)    
    
    # MANUAL
    fig = ana.manual_fitting()
    
    # MINIMIZER
#    mini = ana.minimizer_fitting(check_print=True)
#    ana.model.parameters_0 = mini.x
    
#    # MCMC
#    samp = ana.mcmc_fitting()
    



    