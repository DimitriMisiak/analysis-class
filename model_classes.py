#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: misiak

"""
import numpy as np
import matplotlib.pyplot as plt

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
        
        return fig


if __name__ == '__main__':
    
    plt.close('all')
    plt.rcParams['text.usetex']=True
    
    model = Model_ruo2(model='mmr3')

    model.expo_plot()
