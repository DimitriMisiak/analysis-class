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
import matplotlib.patheffects as pe

from lighten_color import lighten_color

def plot_thresh_cut(xy_data, xy_labels, cut_array, thresh, num):

    x_data, y_data = xy_data
    xlabel, ylabel = xy_labels
    
    x_line = np.linspace(x_data.min(),
                         x_data.max(),
                         100)
    
    fig = plt.figure(num)
    ax = fig.subplots()
    
    ax.loglog(x_data, y_data,
              label='all trig',
              ls='none', marker='+', color='r', alpha=0.1)

    ax.loglog(x_data[cut_array],
              y_data[cut_array],
              label='trig passing cut',
              ls='none', marker='+', color='slateblue')

    if thresh is not None:
        thresh_line = thresh * np.ones(x_line.shape)
        ax.plot(x_line, thresh_line,
                label='Threshold {}={}'.format(ylabel, thresh),
                color='k')
    
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(num)
    ax.legend()
    ax.grid(True)
    
    return fig, ax


def custom_bin_edges(bmin, bmax, res):
    num = int(abs(bmax-bmin)/res)
    return np.linspace(bmin, bmax, num)


def ax_hist(axis, bin_edges, data_array, lab, color='slateblue'):
    """ Draw pretty histogramm and cdf in given axis.
    Return bin_array, hist_array, cdf_array.    
    """
    c_dark = lighten_color(color, 1.5)
    c_light = lighten_color(color, 0.8)
    
    style = [pe.Normal(), pe.withStroke(foreground='k', linewidth=3)]

    bin_array = bin_edges[1:] - (bin_edges[1]-bin_edges[0])/2

    data_hist, _ = np.histogram(data_array, bins=bin_edges)
    
    data_sorted = np.sort(data_array)
    ndim = data_sorted.shape[-1]
    cdf = (np.arange(ndim)+1) / float(ndim)    
    
    hist_line, = axis.plot(bin_array, data_hist, ls='steps-mid',
                           color=c_dark)

    axis.fill_between(bin_array, data_hist, label=lab,
                      color=color, step='mid')

    a0 = axis.twinx()
    a0.set_ylabel('CDF', color='grey')
    a0.tick_params(axis='y', labelcolor='grey')
    
    cdf_line, = a0.plot(data_sorted, cdf,
                        ls='steps', color=c_light, path_effects=style)
    
    axis.grid(True)
    axis.set_ylabel('Counts Events {}'.format(lab), color='k')
    axis.tick_params(axis='y', labelcolor='k')
    axis.set_xlabel('Energy [ADU]')
    
    axis.legend()
    
    axis.set_yscale('log')
    axis.set_xlim(bin_edges[0], bin_edges[-1])
    
    return bin_array, data_hist, cdf


class Artifact(object):
    """ Empty class. To define namespace.
    """
    pass


class Copse(object):
    """ Defining the cuts in a kinda clean/handy way.
    """

    def new_cut(self, name, thruth_array):            
        setattr(self, name, thruth_array)
        

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

        # definng quality cuts
        self.all.trig.cut = Copse()
        self.all.noise.cut = Copse()


#%%
if __name__ == '__main__':
    
    plt.close('all')
    
    run = 'tg10l005'
    ana = Analysis(run)

    trig = ana.all.trig
    noise = ana.all.noise
    
    ### CUT etype EVENTS
    thresh_chi2_heat = 300
    thresh_chi2_ion = 300
    thresh_offset_ion = 14000
    
    etypes = (noise, trig)
    etype_labels = ('noise', 'trig')
    
    for etype, elab in zip(etypes, etype_labels):
    
        energy = etype.filt_decor.Energy_OF
        
        chi2 = etype.filt_decor.chi2_OF
        offset = etype.raw.Off
    
        # CUT Chi2 heat
        etype.cut.new_cut('chi2_heat', chi2[:, 0]<thresh_chi2_heat)
        
        # CUT Chi2 Ion
        etype.cut.new_cut('chi2_ion', np.all(chi2[:, 2:]<thresh_chi2_ion, axis=1))    
        
        # CUT Offset Ion
        etype.cut.new_cut('offset_ion', np.all(offset[:, 2:]<thresh_offset_ion, axis=1))      
        
        # CUT Quality (all cuts)
        quality_truth_array = np.all((etype.cut.chi2_heat,
                                   etype.cut.chi2_ion,
                                   etype.cut.offset_ion), axis=0)
        etype.cut.new_cut('quality', quality_truth_array)
        etype.nsamples_quality = np.count_nonzero(etype.cut.quality)
        
        # chi2 vs Energy plot
        ax_titles = ('Heat', 'Ion A', 'Ion B', 'Ion C', 'Ion D')       
        ax_tuples = ((1, 0), (0, 1), (0, 2), (1, 1), (1, 2))       
        data_ind = (0, 2, 3, 4, 5)       
        x_datas = (np.abs(energy[:, i]) for i in data_ind)    
        y_datas = (chi2[:, i] for i in data_ind)
    
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(11.69, 8.27),
                                 num='{} : Quality Cut Plot'.format(elab))

        for tupl, xdata, ydata, title in zip(ax_tuples, x_datas, y_datas, ax_titles):
            
            ax = axes[tupl]
            
            ax.plot(xdata, ydata,
                    label='All events: {}'.format(etype.nsamples),
                    c='red', marker=',', ls='none')
            
            xdata_cut = xdata[etype.cut.quality]
            ydata_cut = ydata[etype.cut.quality]
            
            ax.plot(xdata_cut, ydata_cut,
                    label='Quality events: {}'.format(etype.nsamples_quality),
                    c='slateblue', marker=',', ls='none')
        
            ax.legend()
            ax.set_title(title)
            ax.set_xlabel('Energy [ADU]')
            ax.set_ylabel('$\chi^2$')
            ax.set_yscale('log')
            ax.set_xscale('log')
            
        fig.delaxes(axes[0,0])    
        fig.tight_layout()


        # Histogramm
        ax_titles = ('Heat', 'Ion A', 'Ion B', 'Ion C', 'Ion D')       
        ax_tuples = ((0, 1), (1, 0), (1, 1), (2, 0), (2, 1))
        data_ind = (0, 2, 3, 4, 5)  
        x_datas = (energy[:, i] for i in data_ind)    
    
        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(11.69, 8.27),
                                 num='{} : Quality Cut Histogram'.format(elab))
        
        for tupl, xdata, title in zip(ax_tuples, x_datas, ax_titles):
            
            ax = axes[tupl]
            
            bin_edges = np.histogram_bin_edges(xdata[etype.cut.quality], bins=250)
        
            ax_hist(ax, bin_edges, xdata,
                    'All events', color='coral')
            ax_hist(ax, bin_edges, xdata[etype.cut.quality],
                    'Quality events', color='slateblue')
            
            ax.set_title(title)
    
        fig.delaxes(axes[0,0])    
        fig.tight_layout()

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
