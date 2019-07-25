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
import abc

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

class Artifact(object):
    """ Empty class. To define namespace.
    """
    pass

class Tree(object):
    """ Tree class, mimicking ROOT's Tree.
    """
    
    def __init__(self, root, tree_name):
        
        self._tree = root[tree_name]
        
        self._keys = list()
        for k in self._tree.keys():
            # sanitize keys and values
            key = k.decode('utf-8')
            (self._keys).append(key)
            
            value = np.array(self._tree[k].array())
            if value.ndim > 1 and value.shape[0]==1:
                value = value[0]
            setattr(self, key, value)

class Wisp(object):
    """ Preparing inheritance for Root and Thicket classes.
    """
    __metaclass__ = abc.ABCMeta
    
    def __init__(self):

        self.run_tree = Artifact()
        
        self._event_types = ('trig', 'noise')
        self._processing_types = ('raw', 'filt', 'filt_decor')   

        for etype in self._event_types:
            setattr(self, etype, Artifact())
            
            for ptype in self._processing_types:
                arti = getattr(self, etype)
                setattr(arti, ptype, Artifact())
        

class Root(Wisp):
    """ Root class, mimicking ROOT's Root.
    """
    def __init__(self, file_path):
        
        Wisp.__init__(self)

        self._root = uproot.open(file_path)
        
        # general info about the run
        self.run_tree = Tree(self._root, 'RunTree_Normal')
        
        for etype in self._event_types:
            arti = getattr(self, etype)
            
            for ptype in self._processing_types:
                tree = Tree(self._root,
                            'EventTree_{}_Normal_{}'.format(etype, ptype))
                setattr(arti, ptype, tree)


class Guardian(Wisp):
    """ Concanetation of Root partitions.
    """
    def __init__(self, root_list):
        
        Wisp.__init__(self)
        
        self.roots = root_list
        
        # self.run_tree for
        # stacking all the run_tree attributes
        # by listing then numpy.stack-ing
        # and
        # self.data_types.tree_labels for
        # concatenate all the data_types.tree_labels attributes
        # by listing then numpy.concatenate-ing
        
        # initializing these attributes with empty lists
        root0 = self.roots[0]
        self.run_tree._keys = root0.run_tree._keys
        for key in self.run_tree._keys:
            setattr(self.run_tree, key, list()) 

        for etype in self._event_types:
            root_arti = getattr(root0, etype)
            guard_arti = getattr(self, etype)
            
            for ptype in self._processing_types:
                root_tree = getattr(root_arti, ptype)
                guard_tree = getattr(guard_arti, ptype)
                setattr(guard_tree, '_keys', root_tree._keys)

                for key in guard_tree._keys:                
                    setattr(guard_tree, key, list())            
        
        # appending the attributes of all root files to the empty list
        for root in self.roots:
            # for self.run_tree
            for key in self.run_tree._keys:
                root_attr = getattr(root.run_tree, key)
                attr_list = getattr(self.run_tree, key)
                attr_list.append(root_attr)

            for etype in self._event_types:
                root_arti = getattr(root, etype)
                guard_arti = getattr(self, etype)
                
                for ptype in self._processing_types:
                    root_tree = getattr(root_arti, ptype)
                    guard_tree = getattr(guard_arti, ptype)
                    
                    for key in guard_tree._keys:
                        root_attr = getattr(root_tree, key)
                        guard_attr = getattr(guard_tree, key)                 
                        guard_attr.append(root_attr)                      
                
        # stacking for self.run_tree
        for key in self.run_tree._keys:
            attr_list = getattr(self.run_tree, key)
            attr_array = np.stack(attr_list, axis=0)
            setattr(self.run_tree, key, attr_array)
   
        # concatenating for self.data_types.tree_labels
        for etype in self._event_types:
            arti = getattr(self, etype)
            
            for ptype in self._processing_types:
                tree = getattr(arti, ptype)
                    
                for key in tree._keys:
                    attr_list = getattr(tree, key)
                    attr_array = np.concatenate(attr_list, axis=0)
                    setattr(tree, key, attr_array)
        
            arti.nsamples = attr_array.shape[0]
                
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
        
        self._part_label = ['part_'+ num for num in self._part_num]
        
        self.roots = list()
        
        for fname, label in zip(self.files, self._part_label):
            fpath = '/'.join((self.dir, fname))
            part_root = Root(fpath)
#            setattr(self, label, part_root)
            (self.roots).append(part_root)

        self.all = Guardian(self.roots)

        self.all.trig.ref = self.all.trig.filt_decor
        

#        amp_list = list()
#        chi2_list = list()
#        amp_noise_list = list()
#        chi2_noise_list = list()
#        for lab in self.part_label:
#            root = getattr(self, lab)
#            amp_list.append(root.trig.filt_decor.Energy_OF)
#            chi2_list.append(root.trig.filt_decor.chi2_OF)
#            amp_noise_list.append(root.noise.filt_decor.Energy_OF_t0)
#            chi2_noise_list.append(root.noise.filt_decor.chi2_OF_t0)
#        
#        self.amp = np.concatenate(amp_list, axis=0)
#        self.chi2 = np.concatenate(chi2_list, axis=0)
#        self.trig_samples = (self.amp).shape[0]
#
#        self.amp_noise = np.concatenate(amp_noise_list, axis=0)
#        self.chi2_noise = np.concatenate(chi2_noise_list, axis=0)
#        self.noise_samples = (self.amp_noise).shape[0]

        self._thresh_chi2_heat = np.inf
        self._cut_chi2_heat_trig = np.ones(self.all.trig.nsamples, dtype=bool)
        self._cut_chi2_heat_noise = np.ones(self.all.noise.nsamples, dtype=bool)

        self._thresh_chi2_ion = np.inf
        self._cut_chi2_ion_trig = np.ones(self.all.trig.nsamples, dtype=bool)
        self._cut_chi2_ion_noise = np.ones(self.all.noise.nsamples, dtype=bool)
#        
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
        
        condi_trig = self.all.trig.filt_decor.chi2_OF[:, 0] < self.thresh_chi2_heat
#        self._cut_chi2_heat_trig = np.nonzero(condi_trig)[0]
        self._cut_chi2_heat_trig = condi_trig
        
        condi_noise = self.all.noise.filt_decor.chi2_OF[:, 0] < self.thresh_chi2_heat
#        self._cut_chi2_heat_noise = np.nonzero(condi_noise)[0]
        self._cut_chi2_heat_noise = condi_noise
 
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
        
        condi_trig = self.all.trig.filt_decor.chi2_OF[:, 2:] < self.thresh_chi2_ion
#        self._cut_chi2_ion_trig = np.nonzero(condi_trig)[0]
        self._cut_chi2_ion_trig = condi_trig
        
        condi_noise = self.all.noise.filt_decor.chi2_OF[:, 2:] < self.thresh_chi2_ion
#        self._cut_chi2_ion_noise = np.nonzero(condi_noise)[0]
        self._cut_chi2_ion_noise = condi_noise
#    

if __name__ == '__main__':
    
    plt.close('all')
    
    fp = '/home/misiak/Data/data_run57/tg10l005/RED80/ProcessedData_tg10l005_S05_RED80_ChanTrig0.root'
    fp2 = '/home/misiak/Data/data_run57/tg10l005/RED80/ProcessedData_tg10l005_S04_RED80_ChanTrig0.root'
        
#    root = Root(fp)
#    print(root.trig.filt_decor.Energy_OF.shape)
#    
#    root2 = Root(fp2)
#    print(root2.trig.filt_decor.Energy_OF.shape)
#    
    run = 'tg10l005'
    ana = Analysis(run)
    
#    gua = Guardian(ana.roots)
    
    # CUT Chi2 heat
    ana.set_cut_chi2_heat(300)
    
    fig, ax = plot_thresh_cut(
            (ana.all.trig.filt_decor.Energy_OF[:, 0], ana.all.trig.filt_decor.chi2_OF[:, 0]),
            ('Amplitude [ADU]', '$\chi_2$'),
            ana.cut_chi2_heat_trig,
            ana.thresh_chi2_heat,
            'Heat Channel: Chi2(amp)'
    )
    
    # CUT Chi2 Ion
    
    ana.set_cut_chi2_ion(300)
    
    fig, ax = plot_thresh_cut(
            (ana.all.trig.filt_decor.Energy_OF[:, 2], ana.all.trig.filt_decor.chi2_OF[:, 2]),
            ('Amplitude [ADU]', '$\chi_2$'),
            ana.cut_chi2_ion_trig[:, 0],
            ana.thresh_chi2_ion,
            'Ion Channel: Chi2(amp)'
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
