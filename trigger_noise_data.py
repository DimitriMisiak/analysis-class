#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 10:30:13 2019

@author: misiak
"""

import numpy as np
import matplotlib.pyplot as plt
import uproot

from core_classes import Tree

###
plt.close('all')

#file_path = '/home/misiak/Data/data_run57/tg21l000/RED70/RootFiles/TriggerData_tg21l000_S02_RED70_ChanTrig0.root'
file_path = '/home/misiak/Data/data_run57/tg21l000/RED70/RootFiles/NoiseData_tg21l000_S02_RED70_ChanTrig0.root'

root = uproot.open(file_path)

tree = Tree(root, 'tree')

tree.samples = tree.Event_Number.shape[0]

#we need some cut
cut = np.all(np.abs(tree.Trace_Heat_A_Raw) < 10000, axis=1)

trace_heat = tree.Trace_Heat_A_Filt_Decor.T[:, cut]
trace_ion = tree.Trace_Ion_B_Raw.T[:, cut]

#Trace Heat Raw
plt.figure('Heat Raw')
plt.plot(trace_heat, alpha=0.1)

##Trace Ion A Raw
#plt.figure('Ion A Raw')
#plt.plot(trace_ion, alpha=0.1)

##Trace OF, meh this is the OF filter
#plt.figure('trace OF')
#plt.plot(tree.Trace_OF.T, alpha=0.1)