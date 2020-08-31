#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: misiak

"""
import numpy as np
import matplotlib.pyplot as plt

<<<<<<< HEAD
import ethem as eth
from data_classes import Data_vi, plot_vi, chi2_list, noising_data

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
                self.system, self.all_parameters, self.eval_dict
        )
    
        self.parameters_0 = [self.eval_dict[p] for p in self.parameters]
    
=======
from ruo2_equation import mmr3_eq_resistance, hart_eq_resistance


class Model_ruo2(object):
#    
#    def function(self, param, temp_array):
#        a1, a2, a3 = param
#        res_array = a2 * np.exp( (a1/temp_array)**(1./a3) )
#        return res_array    

    def _function_mmr3(self, param, temp_array):
        return mmr3_eq_resistance(param, temp_array)
    
    def _function_hart(self, param, temp_array):
        return hart_eq_resistance(param, temp_array)
>>>>>>> f7bb10be8143c9b22931c3c9885b7aac22bae186
    
    def __init__(self, model='mmr3'):

        if model.lower() == 'mmr3':
            self.function = self._function_mmr3
        elif model.lower() == 'hart':
            self.function = self._function_hart
        else:
            raise Exception((
                    'The model \"{}\" is not recognized. '
                    'Available: mmr3, hart'
            ).format(model))
            
<<<<<<< HEAD
            model_iv = list()
            for curr in i_array:
                
                p_top = (temp, curr)
                p_full = p_top + tuple(param)
                
                model_iv.append(self.steady_state_fun(p_full).x[-1])
                
            model_iv_array = np.vstack((i_array, model_iv))
            model_iv_list.append(model_iv_array)
    
        return model_iv_list   
        

    def fake_data(self, param, sigma_coeff=0.1):
=======
        self.parameters_0 = [0.7, 1e+03, 3]
        self.temp_array_0 = np.linspace(12e-3, 100e-3, 100)
        self.res_array_0 = self.function(self.parameters_0, self.temp_array_0)
        self.std_array_0 = 0.1*self.res_array_0
        self.fake_array_0 = self.fake_data(
                self.parameters_0,
                self.temp_array_0,
                np.random.normal(loc=0, scale=self.std_array_0)
        )


    def fake_data(self, param, temp_array, noise_array):
>>>>>>> f7bb10be8143c9b22931c3c9885b7aac22bae186
        """ Add a gaussian noise depending on the model_data.
        """
        model_array = self.function(param, temp_array)
        fake_array = model_array + noise_array
        return fake_array


    def expo_plot(self, num='Model ruo2 expo plot'):
        fig = plt.figure(num=num)
        ax = fig.subplots()
        ax.set_title(num)
        ax.plot(
                self.temp_array_0, 
                self.res_array_0,
                color='slateblue',
                label='model parameters0\n{}'.format(self.parameters_0)
        )
        ax.errorbar(
                self.temp_array_0,
                self.fake_array_0,
                yerr=self.std_array_0,
                ls='none',
                marker='.',
                color='k',
                label='fake data\n(0.1 relative error)'
        )
        ax.set_xlabel('Temperature [K]')
        ax.set_ylabel('Resistance [$\Omega$]')
        ax.legend()
        ax.grid(True)
        fig.tight_layout()
        
<<<<<<< HEAD
        model_data = self.function(param)
        
        fake_data = noising_data(model_data, sigma_coeff)
    
        return fake_data


    def comparator(self, param,
                   data_arrays=None, data_errors=None, 
                   check_print=False):
        
        if data_arrays is None:
            data_arrays = self.data.vi_arrays
        if data_errors is None:
            data_errors = self.data.std_arrays
        model_arrays = self.function(param)
        
        x2 = chi2_list(data_arrays, model_arrays, data_errors)

        if check_print is True:
            print(x2)

        return x2
    

=======
        return fig


>>>>>>> f7bb10be8143c9b22931c3c9885b7aac22bae186
if __name__ == '__main__':
    
    plt.close('all')
    plt.rcParams['text.usetex']=True
    
<<<<<<< HEAD
    cryo = syst.Thermostat_b
    capa = syst.Capacitor_f
    elntd = syst.Resistor_ntd
    leak = syst.ThermalLink_leak
    epcoup = syst.ThermalLink_ep
    glue = syst.ThermalLink_glue

    # "primary" parameters
    model.conditionnal_parameters = (cryo.temperature, capa.current,)

#    # "secondary" parameters
    model.parameters = ()
    
    model.parameters = (
            elntd.R0,
            elntd.T0,
            leak.cond_alpha,
            epcoup.cond_alpha,
    )

    model.lock_model()
    
    modvi_arrays = model.function(model.parameters_0)
    
    fig = plot_vi(model.temp_list, modvi_arrays, None, num=model.data.label)
    
    fake_arrays = model.fake_data(model.parameters_0, 0.1)
    
    plot_vi(model.temp_list, fake_arrays, 0.1, num=model.data.label,
            ls='none', marker='.')
    
    chi2 = chi2_list(modvi_arrays, fake_arrays, 0.1)
    
    plt.grid()
    plt.xlabel('Voltage [V]')
    plt.ylabel('Current [A]')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(title='$\chi^2=${:.3e} ; dof={}'.format(chi2, model.data.nsamples),
               loc='upper left', bbox_to_anchor=(1, 1), ncol=2)
    plt.title(fig.get_label())
    plt.tight_layout(rect=(0., 0., 1., 1.))    
    
    
#    model.data.plot(num=model.data.label, ls='none', marker='.')    
    chi2_data = model.comparator(model.parameters_0)
    chi2_data_bis = chi2_list(modvi_arrays, model.data.vi_arrays, model.data.std_arrays)
    assert chi2_data == chi2_data_bis    
=======
    model = Model_ruo2(model='mmr3')

    model.expo_plot()
>>>>>>> f7bb10be8143c9b22931c3c9885b7aac22bae186
