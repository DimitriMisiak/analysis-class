#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: misiak

"""
import numpy as np
import matplotlib.pyplot as plt

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


class Data_mmr3(object):
    
    def __init__(self, data_path_list, version='old'):
        
        self.data_paths = data_path_list

        self.data_dict = load_mmr3_file(self.data_paths)

        self.time = self.data_dict['Time']
        
        if version == 'old':
            self.temperature = self.data_dict['MMR3-156_1_Conv']
        elif version == 'new':
            self.temperature = self.data_dict['RuO2 MC_Conv']
        else:
            raise Exception(
                    'The value of the keyword "version" is not recognized. '
                    'Choose between "old" and "new".')
        
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
    
    data_paths = ('test/MACRT_2019-07-15.csv', 'test/MACRT_2019-07-16.csv')
    data_mmr3 = Data_mmr3(data_paths, version='old')
    
    data_mmr3.expo_plot()
