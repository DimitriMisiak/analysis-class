#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: misiak

"""

import numpy as np
import matplotlib.pyplot as plt
import mcmc_red as mcr

def load_vi_file(fpath):

    raw_data = (
        np.loadtxt(fpath, skiprows=1, unpack=True)
    )

    temp_array, curr_array, volt_array, resi_array = raw_data
    
    temp_list = np.unique(temp_array)
    
    vi_arrays = list()
    for temp in temp_list:
        
        iv_array = raw_data[1:3, temp_array == temp]
        (vi_arrays).append(iv_array)

    return temp_list, vi_arrays

def plot_vi(temp_list, vi_arrays, std_arrays,
            num='plot vi default num', color=None, mode='errorbar', **kwargs):
    
    ncurves = len(temp_list)
    
    if isinstance(color, (tuple, list, np.ndarray)):
        assert len(color) == ncurves
        
    elif (color is None) or isinstance(color, str):
        color = (color,) * ncurves
        
    else:
        raise Exception('Invalid type for "color" keyword.')
    
    fig = plt.figure(num=num, figsize=(8,5))

    for i in range(ncurves):
        temp = temp_list[i]
        i_array, v_array = vi_arrays[i]
        if std_arrays is None:
            std_array = None
        elif isinstance(std_arrays, float):
            std_array = v_array * std_arrays
        else:
            _, std_array = std_arrays[i]
        c = color[i]
        
        if mode == 'errorbar':
            plt.errorbar(i_array, v_array, yerr=std_array,
                         color=c,
                         label='{0:.1f} mK'.format(temp*1e3),
                         **kwargs)
        
        elif mode == 'plot': 
            plt.plot(i_array, v_array, label='{0:.1f} mK'.format(temp*1e3),
                     color = c,
                     **kwargs)        

    return fig


def chi2_list(data_1, data_2, sigma):
    """ Return the summed chi2 between two given set of data.
    """
    chi2_tot = 0

    # tweak
    if isinstance(sigma, float):
        sigma = [sigma * d for d in data_1]

    for d1, d2, sig in zip(data_1, data_2, sigma):

        f1, y1 = d1
        f2, y2 = d2
        sf, sy = sig
        
        chi2_tot += mcr.chi2_simple(y1, y2, sy)

    return chi2_tot


def blob_noise(vi_arrays, sigma_coeff=0.1):
    
    blob_data = list()
    for md in vi_arrays:

        x_data, y_data = md
        sigma = y_data * sigma_coeff
        blob = np.random.normal(0,sigma,sigma.shape)
        blob_data.append(np.vstack((x_data, blob)))

    return blob_data    


def noising_data(data_arrays, sigma_coeff=0.1):
    
    blob_data = blob_noise(data_arrays, sigma_coeff=sigma_coeff)
    
    noisy_data = list()
    for data_array, blob_array in zip(data_arrays, blob_data):
        
        x_array = data_array[0]
        y_array = data_array[1] + blob_array[1]
        
        noisy_array = np.vstack((x_array, y_array))
        noisy_data.append(noisy_array)
        
    return noisy_data

class Data_vi(object):
    
    data_type = 'VI characteristic'
    
    def __init__(self, data_path, label, error=0.1):
        
        self.data_path = data_path
        self.label = ' '.join((self.data_type, label))

        self.temp_list, self.vi_arrays = load_vi_file(self.data_path)

        if isinstance(error, str):
            temp_list, std_arrays = load_vi_file(error)
            assert temp_list == self.temp_list
            self.std_arrays = std_arrays
       
        else:
            error_coef = float(error)
            self.std_arrays = list()
            for i_array, v_array in self.vi_arrays:
                std_array = np.vstack((i_array, error_coef*v_array))
                self.std_arrays.append(std_array)
        
        nsamples = 0
        for vi_array in self.vi_arrays:
            nsamples += len(vi_array[0])
        self.nsamples = nsamples

    def plot(self, num=None, **kwargs):
        
        if num is None:
            num = self.label
        
        return plot_vi(self.temp_list, self.vi_arrays, self.std_arrays,
                       num=num, **kwargs)
        
if __name__ == '__main__':
    
    plt.close('all')
    plt.rcParams['text.usetex']=True
    
    data_path = 'test/vi_run57_red80_late.csv'
    
    data_vi = Data_vi(data_path, 'RED80 Late')
    
    fig = data_vi.plot(marker='.', lw=0.5)

    plt.grid()
    plt.xlabel('Voltage [V]')
    plt.ylabel('Current [A]')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.title(fig.get_label())
    plt.tight_layout(rect=(0., 0., 1., 1.))
    