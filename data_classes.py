#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: misiak

"""
import numpy as np
import matplotlib.pyplot as plt
import mcmc_red as mcr

def load_mmr3_file(fpath_list):

    header_data = (
            np.loadtxt(fpath_list[0], skiprows=0, unpack=True,
                       delimiter=';', dtype=str, max_rows=1)        
    )
    
    data_dict = dict()
    for title in header_data:
        data_dict[title]=[]
    
    for fpath in fpath_list:
        raw_data = (
            np.loadtxt(fpath, skiprows=1, unpack=True, delimiter=';', dtype=str)
        )
    
        for title, raw_array in zip(header_data, raw_data):
            try:
                data_dict[title]= np.concatenate(
                        (data_dict[title], raw_array.astype(float))
                )
            except:
                data_dict[title]= np.concatenate(
                        (data_dict[title], raw_array)
                )
    
    return data_dict


<<<<<<< HEAD

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
=======
class Data_mmr3(object):
>>>>>>> f7bb10be8143c9b22931c3c9885b7aac22bae186
    
    def __init__(self, data_path_list, version='old'):
        
        self.data_paths = data_path_list

        self.data_dict = load_mmr3_file(self.data_paths)

        self.time = self.data_dict['Time']
        
        if version == 'old':
            self.temperature = self.data_dict['MMR3-156_1_Conv']
        elif version == 'new':
            self.temperature = self.data_dict['RuO2 MC_Conv']
        else:
<<<<<<< HEAD
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
=======
            raise Exception(
                    'The value of the keyword "version" is not recognized. '
                    'Choose between "old" and "new".')
>>>>>>> f7bb10be8143c9b22931c3c9885b7aac22bae186
        
        self.nsamples = self.temperature.shape[0]

    def expo_plot(self, num='Data mmr3 expo plot'):
        fig = plt.figure(num)
        ax = fig.subplots()
        ax.set_title(num) 
        ax.plot(
                self.time,
                self.temperature,
                color='slateblue',
                marker='+',
                label='RuO2 Mixing Chamber'
        )
        ax.set_xlabel('Time Unix [s?]')
        ax.set_ylabel('Temperature MC [K]')
        ax.grid(True)
        ax.legend()
        fig.tight_layout()
        
        return fig
        
if __name__ == '__main__':
    
    plt.close('all')
    plt.rcParams['text.usetex']=True
    
<<<<<<< HEAD
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
    
=======
    data_paths = ('test/MACRT_2019-07-15.csv', 'test/MACRT_2019-07-16.csv')
    data_mmr3 = Data_mmr3(data_paths, version='old')
    
    data_mmr3.expo_plot()
>>>>>>> f7bb10be8143c9b22931c3c9885b7aac22bae186
