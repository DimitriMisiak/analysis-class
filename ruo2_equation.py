#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 16:31:56 2019

@author: misiak
"""
import numpy as np

def mmr3_eq_temperature(param, resistance):
    a, b, c = param
    temperature = a / ( np.log(resistance/b) )**c
    return temperature
   
def mmr3_eq_resistance(param, temperature):
    a, b, c = param
    resistance = b * np.exp( (a/temperature)**(1/c) )
    return resistance

p0 = [0.7, 1000, 3]
diff0 = abs(mmr3_eq_temperature(p0, mmr3_eq_resistance(p0, 0.018)) - 0.018)
assert diff0 < 1e-15


def hart_eq_temperature(param, resistance):
    a, b, c = param
    logr = np.log(resistance)
    temperature = (a + b*logr + c*logr**3)**-1
    return temperature

def hart_eq_resistance(param, temperature):
    a, b, c = param
    y = (a - temperature**-1) / (2*c)
    x = ( ( b/(3*c) )**3 + y**2 )**0.5
    resistance = np.exp( (x-y)**(1./3) - (x+y)**(1./3) )
    return resistance

#p1 = [1e-3, 2e-4, 9e-8]
#diff1 = abs(hart_eq_temperature(p0, hart_eq_resistance(p1, 0.018)) - 0.018)
#assert diff1 < 1e-15


#def limited_eq_temperature(param, resistance):
#    a, b, c = param
#    temperature = 1 / ( 1/c + (1/b)*np.log(resistance/a) )
#    return temperature
#
#def limited_eq_resistance(param, temperature):
#    a, b, c = param
#    resistance = a * np.exp( b * (temperature**-1 - c**-1) )
#    return resistance
#
#p2 = [2e6, 0.001, 0.022]
#diff2 = abs(limited_eq_temperature(p0, limited_eq_resistance(p2, 0.018)) - 0.018)
#assert diff2 < 1e-15
