#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 22:22:22 2019

@author: misiak
"""

import numpy as np

def custom_bin_edges(data_array, res):
    """ Delivers an array of bin edges based on the min and max of the given
    data_array and the pseudo-resolution of the histogramm.    
    """
    bmin = np.min(data_array)
    bmax = np.max(data_array)
    num = int(abs(bmax-bmin)/res)
    
    return np.linspace(bmin, bmax, num)

def cdf_calc(data, axis=-1):
    data_sorted = np.sort(data, axis=-1)
    ndim = data_sorted.shape[-1]
    cdf = (np.arange(ndim)+1) / float(ndim)
    return data_sorted, cdf      