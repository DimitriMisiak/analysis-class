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
import mcmc_red as mcr


def chi2_list(data_1, data_2, sigma):
    """ Return the summed chi2 between two given set of data.
    """
    chi2_tot = 0

    for d1, d2, sig in zip(data_1, data_2, sigma):

        f1, y1 = d1
        f2, y2 = d2
        sf, sy = sig
        
        chi2_tot += mcr.chi2_simple(y1, y2, sy)

    return chi2_tot


def chi2_model(param, check_print=False):
    """ Return the chi2 between the experimental processed data and
    the model evaluated for the given parameter set.
    """
    model_data = make_model_data(param)

    x2 = chi2_list(model_data, exp_data, [0.05*e for e in exp_data])

    if check_print is True:
        print(x2)

    return x2


class Comparator_vi(object):
    
    data_type = 'VI characteristic'
    
    def __init__(self, label):
    
        pass
        


if __name__ == '__main__':
    
    compa = Comparator_vi('compa')
    
    
    