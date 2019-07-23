#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
@author: misiak

"""
import uproot
import numpy as np
import re
import os
import matplotlib.pyplot as plt


def plot_thresh_cut(xy_data, xy_labels, cut_array, thresh, num):

    x_data, y_data = xy_data
    xlabel, ylabel = xy_labels
    
    x_line = np.linspace(x_data.min(),
                         x_data.max(),
                         100)
    
    thresh_line = thresh * np.ones(x_line.shape)
    
    fig = plt.figure(num)
    ax = fig.subplots()
    
    ax.loglog(x_data, y_data,
              label='all trig',
              ls='none', marker='+', color='r', alpha=0.1)

    ax.loglog(x_data[cut_array],
              y_data[cut_array],
              label='trig passing cut',
              ls='none', marker='+', color='slateblue')

    ax.plot(x_line, thresh_line,
            label='Threshold {}={}'.format(ylabel, thresh),
            color='k')
    
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(num)
    ax.legend()
    ax.grid(True)
    
    return fig, ax


class Tree(object):
    """ Tree class
    """
    
    def __init__(self, root, tree_name):
        
        self._tree = root[tree_name]
        
        for k in self._tree.keys():
            # sanitize keys and values
            key = k.decode('utf-8')
            value = np.array(self._tree[k].array())
            if value.ndim > 1 and value.shape[0]==1:
                value = value[0]
            setattr(self, key, value)
            
            
class Root(object):
    """ Root class.
    """
    
    def __init__(self, file_path):
        self._root = uproot.open(file_path)
        
        # general info about the run
        self.run_tree = Tree(self._root, 'RunTree_Normal')
        
        
        self.trig = Tree(self._root, 'EventTree_trig_Normal_filt')
        
        # all processed triggering data from NEPAL
        self.trig_raw = Tree(self._root, 'EventTree_trig_Normal_raw')
        self.trig_filt = Tree(self._root, 'EventTree_trig_Normal_filt')
        self.trig_filt_decor = Tree(self._root, 'EventTree_trig_Normal_filt_decor')
        
        # all processed triggering data from NEPAL
        self.noise_raw = Tree(self._root, 'EventTree_noise_Normal_raw')
        self.noise_filt = Tree(self._root, 'EventTree_noise_Normal_filt')
        self.noise_filt_decor = Tree(self._root, 'EventTree_noise_Normal_filt_decor')

        # main processed data from NEPAL
        ref = 'filt_decor'
        self.trig = getattr(self, 'trig_'+ref)
        self.noise = getattr(self, 'noise_'+ref)


class Analysis(object):
    """ Analysis class.
    """
    
    def __init__(self, run, detector='RED80',
                 run_dir='/home/misiak/Data/data_run57'):
        
        self.run = run
        self.detector = detector
        
        self.dir = '/'.join((run_dir, self.run, self.detector))
    
        self._files_in_data_dir = os.listdir(self.dir)
        
        self.files = []
        for fname in os.listdir(self.dir):
            if 'ProcessedData' in fname:
                self.files.append(fname)
        
        (self.files).sort()
        
        self._part_num = [re.findall('S([0-9]+)_', f)[0] for f in self.files]
        
        self.part_label = ['part_'+ num for num in self._part_num]
        
        for fname, label in zip(self.files, self.part_label):
            fpath = '/'.join((self.dir, fname))
            part_root = Root(fpath)
            setattr(self, label, part_root)

        amp_list = list()
        chi2_list = list()
        amp_noise_list = list()
        chi2_noise_list = list()
        for lab in ana.part_label:
            root = getattr(ana, lab)
            amp_list.append(root.trig.Energy_OF)
            chi2_list.append(root.trig.chi2_OF)
            amp_noise_list.append(root.noise.Energy_OF_t0)
            chi2_noise_list.append(root.noise.chi2_OF_t0)
        
        self.amp = np.concatenate(amp_list, axis=0)
        self.chi2 = np.concatenate(chi2_list, axis=0)
        self.trig_samples = (self.amp).shape[0]

        self.amp_noise = np.concatenate(amp_noise_list, axis=0)
        self.chi2_noise = np.concatenate(chi2_noise_list, axis=0)
        self.noise_samples = (self.amp_noise).shape[0]
      
        self._thresh_chi2_heat = np.inf
        self._cut_chi2_heat_trig = np.ones(self.trig_samples, dtype=bool)
        self._cut_chi2_heat_noise = np.ones(self.noise_samples, dtype=bool)

        self._thresh_chi2_ion = np.inf
        self._cut_chi2_ion_trig = np.ones(self.trig_samples, dtype=bool)
        self._cut_chi2_ion_noise = np.ones(self.noise_samples, dtype=bool)
        
    @property
    def thresh_chi2_heat(self):
        return self._thresh_chi2_heat
    
    @property
    def cut_chi2_heat_trig(self):
        return self._cut_chi2_heat_trig

    @property
    def cut_chi2_heat_noise(self):
        return self._cut_chi2_heat_noise
    
    def set_cut_chi2_heat(self, thresh):
        
        self._thresh_chi2_heat = thresh
        
        condi_trig = self.chi2[:,0] < self.thresh_chi2_heat
        self._cut_chi2_heat_trig = np.nonzero(condi_trig)[0]
        
        condi_noise = self.chi2_noise[:,0] < self.thresh_chi2_heat
        self._cut_chi2_heat_noise = np.nonzero(condi_noise)[0]
 
    @property
    def thresh_chi2_ion(self):
        return self._thresh_chi2_ion
    
    @property
    def cut_chi2_ion_trig(self):
        return self._cut_chi2_ion_trig

    @property
    def cut_chi2_ion_noise(self):
        return self._cut_chi2_ion_noise
    
    def set_cut_chi2_ion(self, thresh):
        
        self._thresh_chi2_ion = thresh
        
        condi_trig = self.chi2[:,0] < self.thresh_chi2_ion
        self._cut_chi2_ion_trig = np.nonzero(condi_trig)[0]
        
        condi_noise = self.chi2_noise[:,0] < self.thresh_chi2_ion
        self._cut_chi2_ion_noise = np.nonzero(condi_noise)[0]    
    

if __name__ == '__main__':
    
    plt.close('all')
    
    run = 'tg10l004'
    ana = Analysis(run)
    
    amp_list = list()
    chi2_list = list()
    
    for lab in ana.part_label:
        root = getattr(ana, lab)
        amp_list.append(root.trig.Energy_OF)
        chi2_list.append(root.trig.chi2_OF)
    
    amp_array = np.concatenate(amp_list, axis=0)
    chi2_array = np.concatenate(chi2_list, axis=0)
    
    # CUT Chi2 heat
    ana.set_cut_chi2_heat(300)
    
    fig, ax = plot_thresh_cut(
            (ana.amp[:, 0], ana.chi2[:, 0]),
            ('Amplitude [ADU]', '$\chi_2$'),
            ana.cut_chi2_heat_trig,
            ana.thresh_chi2_heat,
            'Heat Channel: Chi2(amp)'
    )
    
    # CUT Chi2 Ion
    ana.set_cut_chi2_ion(300)
    
    fig, ax = plot_thresh_cut(
            (ana.amp[:, 0], ana.chi2[:, 0]),
            ('Amplitude [ADU]', '$\chi_2$'),
            ana.cut_chi2_heat_trig,
            ana.thresh_chi2_heat,
            'Heat Channel: Chi2(amp)'
    )
#amp = root.EventTree_noise_Normal_filt_decor.Energy_OF
#chi2 = root.EventTree_noise_Normal_filt_decor.chi2_OF

#        self.run = run_tree['run'].array()[0]
#        self.amp_all = tree_filt_decor['Energy_OF'].array()
#        self.chi2_all = tree_filt_decor['chi2_OF'].array()
#        
##        if amp_array.shape[0] == 0:
##            continue
#    
#        # label non decor = nd
#        amp_array_nd = tree_filt['Energy_OF'].array()
#        chi2_array_nd = tree_filt['chi2_OF'].array()
#    
#        amp_of.append(amp_array)
#        chi2_of.append(chi2_array)
#    
#        cut_ind,_ = full_cut_event(num)
#        
#        amp_of_cut.append(amp_array[cut_ind])
#        chi2_of_cut.append(chi2_array[cut_ind])
#        
#        amp_of_cut_nd.append(amp_array_nd[cut_ind])
#        chi2_of_cut_nd.append(chi2_array_nd[cut_ind])
#        
#        # offset with cuts
#        offset = tree_raw['Off'].array()
#        offset_list.append(offset[cut_ind])
#        
#        # slope with cuts
#        slope = tree_raw['Slope_Ion'].array()
#        slope_list.append(slope[cut_ind])
#        
#        # noise with cuts
#        amp_array_noise = tree_noise['Energy_OF_t0'].array()
#        chi2_array_noise = tree_noise['chi2_OF_t0'].array()
#    
#        cut_ind_noise,_ = full_cut_noise(num)
#        
#        amp_of_noise.append(amp_array_noise[cut_ind_noise])
#        chi2_of_noise.append(chi2_array_noise[cut_ind_noise])    
#        
#        # time stamp with cut
#        micro_stp = tree_filt_decor['MicroStp'].array()
#        num_part = tree_filt_decor['NumPart'].array()
#        freq = run_tree['f_max_heat'].array()[0]
#        hour_stp = (micro_stp/3600.) / freq + num_part
#        
#        pulse_stp.append(hour_stp[cut_ind])
#        run_duration = hour_stp[-1]
#        
#        # gain
#        chan_gain = run_tree['Chan_Gain'].array()[0]
#        gain_chal = chan_gain[0]
#        gain_ion = chan_gain[2:]
#        
#        # maintenance
#    #    maint_cycle = run_tree['MaintenanceCycle'].array()[0] *1.05 / 3600.
#        maint_cycle = run_tree['MaintenanceCycle'].array()[0] / 3600.
#        maint_duration = run_tree['MaintenanceDuration'].array()[0] / 3600.        
#
#    
