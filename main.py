#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
@author: misiak

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from lighten_color import lighten_color

from core_classes import Analysis, Artifact

def custom_bin_edges(data_array, res):
    """ Delivers an array of bin edges based on the min and max of the given
    data_array and the pseudo-resolution of the histogramm.    
    """
    bmin = np.min(data_array)
    bmax = np.max(data_array)
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
                      color=c_light, step='mid')

    a0 = axis.twinx()
    a0.set_ylabel('CDF', color='grey')
    a0.tick_params(axis='y', labelcolor='grey')
    
    cdf_line, = a0.plot(data_sorted, cdf,
                        ls='steps', color=color, path_effects=style)
    
    axis.grid(True)
    axis.set_ylabel('Counts Events {}'.format(lab), color='k')
    axis.tick_params(axis='y', labelcolor='k')
    axis.set_xlabel('Energy [ADU]')
    
    axis.legend(loc=2)
    
    axis.set_yscale('log')
    axis.set_xlim(bin_edges[0], bin_edges[-1])
    
    return bin_array, data_hist, cdf

if __name__ == '__main__':
    
    plt.close('all')
    plt.rcParams['text.usetex']=True
    
    run = 'tg17l003'
    ana = Analysis(run)

    run_tree = ana.all.run_tree
    
    run_tree.chan_label = np.array([
            'heat_a', 'heat_b', 'ion_a', 'ion_b', 'ion_c', 'ion_d'
    ])

    run_tree.chan_valid = np.array([
            0, 2, 3, 4, 5
    ])
    
    run_tree.chan_heat = np.array([0, 1])
    run_tree.chan_ion = np.array([2, 3, 4, 5])
    run_tree.chan_collect = np.array([3, 5])
    run_tree.chan_veto = np.array([2, 4])
    
    trig = ana.all.trig
    noise = ana.all.noise
    
    etypes = (noise, trig)
    
# =============================================================================
# CUT etype EVENTS
# =============================================================================
    thresh_chi2_heat = 300
    thresh_chi2_ion = 300
    thresh_offset_ion = 14000
    
    for etype in etypes:

        if etype is noise:
            energy = etype.filt_decor.Energy_OF_t0
            chi2 = etype.filt_decor.chi2_OF_t0        
        
        if etype is trig:
            energy = etype.filt_decor.Energy_OF
            chi2 = etype.filt_decor.chi2_OF
        
        offset = etype.raw.Off
    
        # CUT Chi2 heat
        etype.cut.new_cut('chi2_heat',
                          chi2[:, 0]<thresh_chi2_heat)
        
        # CUT Chi2 Ion
        etype.cut.new_cut('chi2_ion',
                          np.all(chi2[:, 2:]<thresh_chi2_ion, axis=1))    
        
        # CUT Offset Ion
        etype.cut.new_cut('offset_ion',
                          np.all(offset[:, 2:]<thresh_offset_ion, axis=1))      
        
        # CUT Quality (all cuts)
        quality_truth_array = np.all((etype.cut.chi2_heat,
                                   etype.cut.chi2_ion,
                                   etype.cut.offset_ion), axis=0)
        etype.cut.new_cut('quality', quality_truth_array)
        etype.nsamples_quality = np.count_nonzero(etype.cut.quality)
        
        
#%%
# =============================================================================
# TEMPORAL PLOT
# =============================================================================
    # acquisition frequency
    freq_root = np.ravel((run_tree.f_max_heat, run_tree.f_max_ion))
    assert np.all(freq_root == freq_root[0])
    run_tree.freq = freq_root[0]
    
    # time_array in hours
    secondes_array = trig.filt_decor.MicroStp / run_tree.freq
    hours_array =  secondes_array / 3600 + trig.filt_decor.NumPart 
    run_tree.time = hours_array
    
    # cut on data
    cut = np.ones(run_tree.time.shape, dtype=bool) #no cut
    cut = trig.cut.quality
    
    
    # data
    energy = trig.filt_decor.Energy_OF[cut]
    offset = trig.raw.Off[cut]
    slope = trig.raw.Slope_Ion[cut]
    time = run_tree.time[cut]
    
    # Init figure
    num = 'Temporal Run Check'
    fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(8.27, 11.69),
                             sharex=True, num=num)
    
    # heat trig vs time
    ax = axes[0]
    ax.plot(time, energy[:, 0],
            label='heat a',
            ls='none', marker=',')
    ax.set_ylabel('Energy Heat [ADU]')
    ax.set_yscale('symlog')

    # ion trig vs time
    ax = axes[1]
    for ind in run_tree.chan_ion:
        label = (run_tree.chan_label[ind]).replace('_', ' ')
        ax.plot(time, energy[:, ind],
                label=label,
                ls='none', marker=',')
    ax.set_ylabel('Energy Ion [ADU]')
    ax.set_yscale('symlog')
    
    # heat offset vs time
    ax = axes[2]
    ax.plot(time, offset[:, 0],
            label='heat a',
            ls='none', marker=',')
    ax.set_ylabel('Offset Heat [ADU]')
    
    # ion offset vs time
    ax = axes[3]
    for ind in run_tree.chan_ion:
        label = (run_tree.chan_label[ind]).replace('_', ' ')
        ax.plot(time, offset[:, ind],
                label=label,
                ls='none', marker=',')
    ax.set_ylabel('Offset Ion [ADU]')
    
    # ion slope vs time
    ax = axes[4]
    for ind in run_tree.chan_ion:
        label = (run_tree.chan_label[ind]).replace('_', ' ')
        ax.plot(time, slope[:, ind-2],
                label=label,
                ls='none', marker=',')
    ax.set_ylabel('Slope Ion [ADU/s]')
    
    for ax in axes:
        ax.legend(loc=2)
        ax.grid(True, alpha=0.3)
        
        if ax is not axes[-1]:
            # removing the first tick label
            yticks = ax.yaxis.get_major_ticks()
            yticks[0].label1.set_visible(False)

        if ax is axes[-1]:
            ax.set_xlabel('Time [hours]')

    
    
    fig.text(0.5, 0.98, num,
             horizontalalignment='center',
             verticalalignment='center',
             bbox=dict(facecolor='lime', alpha=0.5))
    
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.subplots_adjust(hspace=.0)

#%%
# =============================================================================
# PLOT Chi2 vs Energy
# =============================================================================
    for etype in etypes:

        if etype is noise:
            energy = etype.filt_decor.Energy_OF_t0
            chi2 = etype.filt_decor.chi2_OF_t0        
        
        if etype is trig:
            energy = etype.filt_decor.Energy_OF
            chi2 = etype.filt_decor.chi2_OF
            
        # chi2 vs Energy plot
        ax_titles = run_tree.chan_label
        ax_tuples = ((1, 0), (0, 1), (0, 2), (1, 1), (1, 2))       
        data_ind = run_tree.chan_valid     
        x_datas = (np.abs(energy[:, i]) for i in data_ind)    
        y_datas = (chi2[:, i] for i in data_ind)
    
        num = '{} : Quality Cut Plot'.format(etype.name)
        
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(11.69, 8.27),
                                 num=num)

        for tupl, xdata, ydata, title in zip(ax_tuples, x_datas, y_datas, ax_titles):
            
            ax = axes[tupl]
            
            ax.plot(xdata, ydata,
                    label='All events: {}'.format(etype.nsamples),
                    c='red', marker=',', ls='none')
            
            xdata_cut = xdata[etype.cut.quality]
            ydata_cut = ydata[etype.cut.quality]
            
            if etype.nsamples < 1000:
                marker = '.'
            else:
                marker = ','

            ax.plot(xdata_cut, ydata_cut,
                    label='Quality events: {}'.format(etype.nsamples_quality),
                    c='slateblue', marker=marker, ls='none')
        
            ax.legend()
            ax.set_title(title.replace('_', ' '))
            ax.set_xlabel('Energy [ADU]')
            ax.set_ylabel('$\chi^2$')
            ax.set_yscale('log')
            ax.set_xscale('log')
            
            ax.set_xlim(xdata_cut.min()*0.5, ax.get_xlim()[1])
            ax.set_ylim(ydata_cut.min()*0.5, ax.get_ylim()[1])

        fig.text(0.5, 0.98, num,
                 horizontalalignment='center',
                 verticalalignment='center',
                 bbox=dict(facecolor='lime', alpha=0.5))
            
        fig.delaxes(axes[0,0])    
        fig.tight_layout(rect=(0, 0, 1, 0.98))

# =============================================================================
# BASELINE RESOLUTION
# =============================================================================
    ana.sigma0 = Artifact('sigma0')
    
    energy = noise.filt_decor.Energy_OF_t0
    
    # for all the channels 
    for ind in run_tree.chan_valid:
        chan = run_tree.chan_label[ind]
        sigma = np.std(energy[:, ind])
        setattr(ana.sigma0, chan, sigma)
    
    # XXX not really good to do the sum just now
    # beacause the energy is still in ADU
    # we need to do the sum with the energy in eV
    
    # for the virtual all_ion channel
    energy_ion_all = np.sum(energy[:, run_tree.chan_ion], axis=0)
    ana.sigma0.ion_all = np.std(energy_ion_all)
    
    # creating the sign correction used when "adding ion channels"
    proto_sign = np.concatenate((run_tree.Sign_Chal, run_tree.Polar_Ion), axis=1)
    assert np.all(proto_sign == proto_sign[0])
    proto_sign = proto_sign[0]
    
    run_tree.sign_corr = np.where(proto_sign<0, -1, 1)
    energy_corr = energy * run_tree.sign_corr
    
    # for the virtual collect channel
    energy_collect = np.sum(energy_corr[:, run_tree.chan_collect], axis=0)
    ana.sigma0.collect = np.std(energy_collect)
    
    # for the virtual veto channel
    energy_veto = np.sum(energy_corr[:, run_tree.chan_veto], axis=0)
    ana.sigma0.veto = np.std(energy_veto)
    
# =============================================================================
# HISTOGRAM
# =============================================================================
    for etype in etypes:

        if etype is noise:
            energy = etype.filt_decor.Energy_OF_t0
        
        if etype is trig:
            energy = etype.filt_decor.Energy_OF
        
        ax_tuples = ((0, 1), (1, 0), (1, 1), (2, 0), (2, 1))
        data_ind = run_tree.chan_valid
    
        num = '{} : Quality Cut Histogram'.format(etype.name)
    
        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(11.69, 8.27),
                                 num=num)
        
        for tupl, ind in zip(ax_tuples, data_ind):
            
            xdata = energy[:, ind]
            label = run_tree.chan_label[ind]
            
            ax = axes[tupl]
            xdata_qual = xdata[etype.cut.quality]
            
            if etype is noise:
                bin_edges = np.histogram_bin_edges(xdata[etype.cut.quality])
            
            if etype is trig:
                bin_edges = custom_bin_edges(xdata_qual, 
                                             getattr(ana.sigma0, label))
        
            ax_hist(ax, bin_edges, xdata,
                    'All events', color='coral')
            ax_hist(ax, bin_edges, xdata_qual,
                    'Quality events', color='slateblue')
            
            ax.set_title(label.replace('_', ' '))
            
            if etype is noise:
                ax.set_yscale('linear')
        
        fig.text(0.5, 0.98, num,
                 horizontalalignment='center',
                 verticalalignment='center',
                 bbox=dict(facecolor='lime', alpha=0.5))
    
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
