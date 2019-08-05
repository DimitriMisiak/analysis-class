#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
@author: misiak

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

import scipy.stats as st

from lighten_color import lighten_color

from core_classes import Analysis, Artifact


import matplotlib.text as mtext
class LegendTitle(object):
    def __init__(self, text_props=None):
        self.text_props = text_props or {}
        super(LegendTitle, self).__init__()

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        title = mtext.Text(x0, y0, r'\underline{' + orig_handle + '}',
                           usetex=True, **self.text_props)
        handlebox.add_artist(title)
        return title
    
from matplotlib.legend import Legend
Legend.update_default_handler_map({str: LegendTitle()})


def custom_bin_edges(data_array, res):
    """ Delivers an array of bin edges based on the min and max of the given
    data_array and the pseudo-resolution of the histogramm.    
    """
    bmin = np.min(data_array)
    bmax = np.max(data_array)
    num = int(abs(bmax-bmin)/res)
    
    return np.linspace(bmin, bmax, num)


def custom_autoscale(axis, xdata, ydata):
    
    gold2 = ((1+5**0.5)/2)**0.5
    
    xmin, xmax = np.min(xdata), np.max(xdata)
    ymin, ymax = np.min(ydata), np.max(ydata)
    
    xcore = (xmax + xmin)/2
    ycore = (ymax + ymin)/2
    
    xside = (xmax - xmin) * gold2
    yside = (ymax - ymin) * gold2
    
    xinf, xsup = xcore - xside/2, xcore + xside/2
    yinf, ysup = ycore - yside/2, ycore + yside/2
    
    axis.set_xlim(xinf, xsup)
    axis.set_ylim(yinf, ysup)


def cfd(data, axis=-1):
    data_sorted = np.sort(data, axis=-1)
    ndim = data_sorted.shape[-1]
    cdf = (np.arange(ndim)+1) / float(ndim)
    return data_sorted, cdf      

def ax_hist(axis, bin_edges, data_array, lab, color='slateblue'):
    """ Draw pretty histogramm and cdf in given axis.
    Return bin_array, hist_array, cdf_array.    
    """
    c_dark = lighten_color(color, 1.5)
    c_light = lighten_color(color, 0.8)
    
    style = [pe.Normal(), pe.withStroke(foreground='k', linewidth=3)]

    bin_array = bin_edges[1:] - (bin_edges[1]-bin_edges[0])/2

    data_hist, _ = np.histogram(data_array, bins=bin_edges)
    
#    data_sorted = np.sort(data_array)
#    ndim = data_sorted.shape[-1]
#    cdf = (np.arange(ndim)+1) / float(ndim)    
    data_sorted, cdf = cfd(data_array)
    
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

class double_norm(st.rv_continuous):
    """ Double Gaussian distribution. """

    def _cdf(self, x, f, loc1, scale1, loc2, scale2):
        cdf1 = (1-f) * st.norm.cdf(x, loc=loc1, scale=scale1)
        cdf2 = f * st.norm.cdf(x, loc=loc2, scale=scale2)
        cdf = cdf1 + cdf2
        return cdf

    def _pdf(self, x, f, loc1, scale1, loc2, scale2):
        pdf1 = (1-f)* st.norm.pdf(x, loc=loc1, scale=scale1)
        pdf2 = f* st.norm.pdf(x, loc=loc2, scale=scale2)
        pdf = pdf1 + pdf2
        return pdf

    def _argcheck(self, *args):
        """Default check for correct values on args and keywords.
        Returns condition array of 1's where arguments are correct and
         0's where they are not.
        """
        cond = 1
#        for arg in args:
#            cond = np.logical_and(cond, (np.asarray(arg) > 0))
        return cond


#%%
if __name__ == '__main__':
    
    # first command
    plt.close('all')
    plt.rcParams['text.usetex']=True
    style = [pe.Normal(), pe.withStroke(foreground='k', linewidth=3)]
    
    run = 'tg25l019'
    ana = Analysis(run)

# =============================================================================
# COMFORT VARIABLES AND ATTRIBUTES
# =============================================================================
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
    
    run_tree.chan_signal = np.append([0,], run_tree.chan_collect)
    
    trig = ana.all.trig
    noise = ana.all.noise
    
    etypes = (noise, trig)

# =============================================================================
# TEMPORAL DATA EXTRACTION
# =============================================================================
    # acquisition frequency
    freq_root = np.ravel((run_tree.f_max_heat, run_tree.f_max_ion))
    assert np.all(freq_root == freq_root[0])
    run_tree.freq = freq_root[0]
    
    # time_array in hours
    def time_stamp(microstp_root, numpart_root):
        secondes_array = microstp_root / run_tree.freq
        hours_array =  secondes_array / 3600 + numpart_root
        return hours_array
    
    # maintenance info
    maint_cycle = run_tree.MaintenanceCycle
    assert np.all(maint_cycle == maint_cycle[0])
    run_tree.maint_cycle = maint_cycle[0] / 3600 # in hours
    
    maint_duration = run_tree.MaintenanceDuration
    assert np.all(maint_duration == maint_duration[0])
    run_tree.maint_duration = maint_duration[0] / 3600 # in hours
    
    # checking if the considered time is within a maintenance
    def maintenance_cut(time_array):
        full_cycle = run_tree.maint_duration + run_tree.maint_cycle
        remainder_array = time_array % full_cycle
        # thruth array, True if time not in a maintenance
        return remainder_array > run_tree.maint_duration

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
        etype.time = time_stamp(etype.filt_decor.MicroStp,
                                etype.filt_decor.NumPart)
    
        # CUT Maintenance time
        etype.cut.new_cut('maintenance',
                          maintenance_cut(etype.time))
    
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
        quality_truth_array = np.all((
                etype.cut.maintenance,
                etype.cut.chi2_heat,
                etype.cut.chi2_ion,
                etype.cut.offset_ion
        ), axis=0)
        etype.cut.new_cut('quality', quality_truth_array)
        
        # number of events passing quality cuts
        etype.nsamples_quality = np.count_nonzero(etype.cut.quality)
        
# =============================================================================
# BASELINE RESOLUTION
# =============================================================================
    noise.sigma0 = Artifact('sigma0')
    
    energy = noise.filt_decor.Energy_OF_t0[noise.cut.quality]
    
    # for all the channels 
    for ind in run_tree.chan_valid:
        chan = run_tree.chan_label[ind]
        sigma = np.std(energy[:, ind])
        setattr(noise.sigma0, chan, sigma)

# =============================================================================
# FIDUCIAL CUT (only on trig)
# =============================================================================
    energy = trig.filt_decor.Energy_OF
    
    cond_veto = []
    
    # fiducial condition
    for ind in run_tree.chan_veto:
        lab = run_tree.chan_label[ind]
        sigma0 = getattr(noise.sigma0, lab)
        
        # consider cut at 2-sigma
        cond = np.abs(energy[:, ind]) < 2*sigma0
    
        cond_veto.append(cond)

    trig.cut.fiducial_raw = np.logical_and(*cond_veto)
    trig.cut.fiducial = np.logical_and(trig.cut.fiducial_raw, trig.cut.quality)
    
    # number of events passing quality and fiducial cut
    trig.nsamples_fiducial = np.count_nonzero(trig.cut.fiducial)

#%%
# =============================================================================
# ESTIMATION OF THE SENSIBILITY
# =============================================================================
    ana.calibration_peak = Artifact('calibration_peak')
    ana.calibration_peak.energy = 10.37 * 1e3 # Ge 10.37keV
    ana.calibration_peak.sigma = Artifact('sigma')
    
    ana.sensitivity = Artifact('sensitivity')
    
    ana.model = Artifact('model')
    ana.model.popt = Artifact('popt')
    
    # check if the run is heat only, and choose a good cut for the data
    if np.all(run_tree.Polar_Ion == 0):
        # heat only run
        ana.calibration_peak.cut_type = 'quality'
    else:
        # heat + ion run
        ana.calibration_peak.cut_type = 'fiducial'
    
    cut_peak = getattr(trig.cut, ana.calibration_peak.cut_type)
    energy = trig.filt_decor.Energy_OF[cut_peak]

    # double gaussian distribution model
    ana.model.dist = double_norm(name='double_gaussian')

    for ind in run_tree.chan_signal:
        
        lab = run_tree.chan_label[ind]
        data = energy[:, ind]
    
        mu01, mu02 = np.quantile(data, [0.10, 0.90])
        sig01 = sig02 = noise.sigma0.heat_a
        p0_light = [0.5, mu01, sig01, mu02, sig02]
        p0 = np.append(p0_light, [0, 1])
    
        popt = ana.model.dist.fit(data, *p0_light, floc=0, fscale=1)
    
        setattr(ana.model.popt, lab, popt)
    
        f, mu1, sig1, mu2, sig2, loc, scale = popt
        gauss1 = (mu1, sig1)
        gauss2 = (mu2, sig2)
        if abs(mu1) > abs(mu2):
            mu, sig = gauss1
        else:
            mu, sig = gauss2
        
        setattr(ana.calibration_peak, lab, mu)
        setattr(ana.calibration_peak.sigma, lab, sig)
        
        sens = getattr(ana.calibration_peak, lab) / ana.calibration_peak.energy
        setattr(ana.sensitivity, lab, sens)

#%%
# =============================================================================
# CONVERSION FROM ADU TO EV
# =============================================================================
    # energy in eV, and energy in eV corrected with the sign
    trig.energy_ev = Artifact('energy_ev')
    noise.energy_ev = Artifact('energy_ev')
    noise.sigma0_ev = Artifact('sigma0_ev')
    ana.calibration_peak.sigma_ev = Artifact('sigma_ev')
    
    # recovering energy in adu
    energy_adu = trig.filt_decor.Energy_OF
    noise_energy_adu = noise.filt_decor.Energy_OF_t0
    sigma0 = noise.sigma0
    sigma = ana.calibration_peak.sigma

    # creating the sign correction used when "adding ion channels"
    proto_sign = np.concatenate((-run_tree.Sign_Chal, run_tree.Polar_Ion), axis=1)
    assert np.all(proto_sign == proto_sign[0])
    proto_sign = proto_sign[0]
    
    # no correction for negative polarisation collecting positive signal,
    # -1 correction for the positive polarisation collecting negative signal
    run_tree.sign_corr = np.where(proto_sign<0, +1, -1)

    for ind in run_tree.chan_signal:
        lab = run_tree.chan_label[ind]
        
        e_adu = energy_adu[:, ind]
        noise_e_adu = noise_energy_adu[:, ind]
        sig0_adu = getattr(sigma0, lab)
        sig_adu = getattr(sigma, lab)
        
        sens = getattr(ana.sensitivity, lab)
        sign = run_tree.sign_corr[ind]        

        e_ev = e_adu / sens
        noise_e_ev = noise_e_adu / sens
        sig0_ev = sign * sig0_adu / sens
        sig_ev = sign * sig_adu / sens

        setattr(trig.energy_ev, lab, e_ev)
        setattr(noise.energy_ev, lab, noise_e_ev)
        setattr(noise.sigma0_ev, lab, sig0_ev)
        setattr(ana.calibration_peak.sigma_ev, lab, sig_ev)

# =============================================================================
# VIRTUAL COLLECT CHANNEL
# =============================================================================
    # labels of the channel with energy conversion in eV
    run_tree.chan_label_virtual = ['heat_a', 'collect']
    
    energy_collect = list()
    noise_collect = list()
    for ind in run_tree.chan_collect:
        lab = run_tree.chan_label[ind]
        run_tree.chan_label_virtual.append(lab)
        
        energy = getattr(trig.energy_ev, lab)
        energy_collect.append(energy)
        
        noise_energy = getattr(noise.energy_ev, lab)
        noise_collect.append(noise_energy)
    
    trig.energy_ev.collect = np.sum(energy_collect, axis=0)
    noise.energy_ev.collect = np.sum(noise_collect, axis=0)

    noise.sigma0_ev.collect = np.std(noise.energy_ev.collect)
## =============================================================================
## DEFINITON OF THE VIRTUAL ELECTRODE (SUMMING THE DATA)
## =============================================================================
#    # XXX not really good to do the sum just now
#    # beacause the energy is still in ADU
#    # we need to do the sum with the energy in eV
#    
#    # for the virtual all_ion channel
#    energy_ion_all = np.sum(energy[:, run_tree.chan_ion], axis=0)
#    noise.sigma0.ion_all = np.std(energy_ion_all)
#    
#    # for the virtual veto channel
#    energy_veto = np.sum(energy_corr[:, run_tree.chan_veto], axis=0)
#    noise.sigma0.veto = np.std(energy_veto)
#        

#%%
# =============================================================================
# TEMPORAL PLOT
# =============================================================================
    
    # cut on data
    cut = np.ones(trig.time.shape, dtype=bool) #no cut
    cut = trig.cut.quality
    
    # data
    energy = trig.filt_decor.Energy_OF
    chi2 = trig.filt_decor.chi2_OF
    offset = trig.raw.Off
    slope = trig.raw.Slope_Ion
    time = trig.time
    
    # Init figure
    num = 'Temporal Run Check'
    fig, axes = plt.subplots(nrows=6, ncols=1, figsize=(8.27, 11.69),
                             sharex=True, num=num)

    # heat trig vs time
    ax = axes[0]
    ax.set_ylabel('Energy Heat [ADU]')
    ax.set_yscale('symlog')
    
    ax.plot(
            time[cut], energy[cut, 0],
            label='heat a', zorder=10,
            ls='none', marker='2', mew=0.8,
    )
    ax.autoscale(False)
    ax.plot(
            time, energy[:, 0],
            label='All events',
            ls='none', marker=',', color='silver',
    )

    # ion trig vs time
    ax = axes[1]
    ax.set_ylabel('Energy Ion [ADU]')
    ax.set_yscale('symlog')
    
    for i, ind in enumerate(run_tree.chan_ion):
        label = (run_tree.chan_label[ind]).replace('_', ' ')
        ax.plot(
                time[cut], energy[cut, ind],
                label=label, zorder=10,
                ls='none', marker=str(i+1), mew=0.8
        )
    ax.autoscale(False)
    ax.plot(
            time, energy[:, run_tree.chan_ion],
            label='All events',
            ls='none', marker=',', color='silver',
    )

    # heat offset vs time
    ax = axes[2]
    ax.plot(
            time[cut], offset[cut, 0],
            label='heat a', zorder=10,
            ls='none', marker='2'
    )
    ax.autoscale(False)
    ax.plot(
            time, offset[:, 0],
            label='All events',
            ls='none', marker=',', c='silver'
    )
    ax.set_ylabel('Offset Heat [ADU]')
    
    # ion offset vs time
    ax = axes[3]
    ax.set_ylabel('Offset Ion [ADU]')
    for i, ind in enumerate(run_tree.chan_ion):
        label = (run_tree.chan_label[ind]).replace('_', ' ')
        ax.plot(
                time[cut], offset[cut, ind],
                label=label, zorder=10,
                ls='none', marker=str(i+1), mew=0.8
        )
    ax.autoscale(False)
    ax.plot(
            time, offset[:, ind],
            label='All events',
            ls='none', marker=',', c='silver'
    )
    
    # ion slope vs time
    ax = axes[4]
    ax.set_ylabel('Slope Ion [ADU/s]')
    for i, ind in enumerate(run_tree.chan_ion):
        label = (run_tree.chan_label[ind]).replace('_', ' ')
        ax.plot(time[cut], slope[cut, ind-2],
                label=label, zorder=10,
                ls='none', marker=str(i+1),
        )
    ax.autoscale(False)
    ax.plot(time, slope[:, ind-2],
            label='All events',
            ls='none', marker=',', c='silver'
    )

    # chi2 vs time
    ax = axes[5]
    ax.set_ylabel('$\chi^2$')
    label = 'chi2 heat A'
    ax.plot(time[cut], chi2[cut, 0],
            label=label, zorder=10,
            ls='none', marker=str(i+1),
    )
    ax.autoscale(False)
    ax.plot(time, chi2[:, 0],
            label='All events',
            ls='none', marker=',', c='silver'
    )
    
    # formatting the axes
    for ax in axes:
        ax.grid(True, alpha=0.3)
        
        # custom legend
        handles = ['Quality events:',]
        labels = ['',]
        for line in ax.get_lines():
            label = line.get_label()
            if label == 'All events':
                if label != labels[0]:
                    handles.insert(0, line)
                    labels.insert(0, label)
            else:
                handles.append(line)
                labels.append(label)
            
        ax.legend(
                handles, labels, loc=2, framealpha=1,
                bbox_to_anchor=(1.05, 1), borderaxespad=0.,
        )

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

#%%
# =============================================================================
# HISTOGRAM ADU
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
                                             getattr(noise.sigma0, label))
        
            ax_hist(ax, bin_edges, xdata,
                    'All events', color='coral')
            ax_hist(ax, bin_edges, xdata_qual,
                    'Quality events', color='slateblue')
            
            if etype is trig:
                xdata_fid = xdata[trig.cut.fiducial]
                ax_hist(ax, bin_edges, xdata_fid,
                        'Fiducial events', color='limegreen')                
                
                if ind in run_tree.chan_signal:
                    popt = getattr(ana.model.popt, label)
                    xrange = np.linspace(xdata_fid.min(), xdata_fid.max(), 1000)
                    pdf = ana.model.dist.pdf(xrange, *popt)
                    
                    normalization = getattr(trig,
                                            'nsamples_{}'.format(
                                                    ana.calibration_peak.cut_type
                                            ))
                    pdf_norm = pdf * normalization * (bin_edges[1] - bin_edges[0])
                    
                    ax.autoscale(False)
                    ax.plot(xrange, pdf_norm,
                            ls='--', color='yellow',
                            label='fit')
            
            ax.legend(loc=2)
            ax.set_title(label.replace('_', ' '))
            
            if etype is noise:
                ax.set_yscale('linear')
        
        fig.text(0.5, 0.98, num,
                 horizontalalignment='center',
                 verticalalignment='center',
                 bbox=dict(facecolor='lime', alpha=0.5))
    
        fig.delaxes(axes[0,0])    
        fig.tight_layout()

#%%
# =============================================================================
# HISTOGRAM EV
# =============================================================================
    for etype in etypes:

        energy = etype.energy_ev
        
        ax_tuples = ((0, 0), (0, 1), (1, 0), (1, 1))
        labels = run_tree.chan_label_virtual
    
        num = '{} : Quality Cut Histogram EV'.format(etype.name)
    
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(11.69, 8.27),
                                 num=num)
        
        for tupl, label in zip(ax_tuples, labels):
            xdata = getattr(energy, label)
            ax = axes[tupl]
            xdata_qual = xdata[etype.cut.quality]
            
            if etype is noise:
                bin_edges = np.histogram_bin_edges(xdata[etype.cut.quality])
            
            if etype is trig:
                bin_edges = custom_bin_edges(xdata_qual, 
                                             getattr(noise.sigma0_ev, label))
        
            ax_hist(ax, bin_edges, xdata,
                    'All events', color='coral')
            ax_hist(ax, bin_edges, xdata_qual,
                    'Quality events', color='slateblue')
            
            if etype is trig:
                xdata_fid = xdata[trig.cut.fiducial]
                ax_hist(ax, bin_edges, xdata_fid,
                        'Fiducial events', color='limegreen')     
                
            ax.set_xlabel('Enregy [eV]')
            ax.legend(loc=2)
            ax.set_title(label.replace('_', ' '))
            
            if etype is noise:
                ax.set_yscale('linear')
        
        fig.text(0.5, 0.98, num,
                 horizontalalignment='center',
                 verticalalignment='center',
                 bbox=dict(facecolor='lime', alpha=0.5))
      
        fig.tight_layout()

#%%
# =============================================================================
# ION VS ION
# =============================================================================

    # recovering data
    energy = trig.filt_decor.Energy_OF
    cut_qual = trig.cut.quality
    cut_fid = trig.cut.fiducial
    
    # initializing pseudo-corner plot
    ax_tuples = [(0,0), (1,0), (1,1), (2,0), (2,1), (2,2)]
    ax_discard = [(0, 1), (1, 2), (0, 2)]
    
    chan_x = np.insert(run_tree.chan_veto, 0, run_tree.chan_collect[1])
    chan_y = np.append(run_tree.chan_veto, run_tree.chan_collect[0])
    
    num = 'Ion vs Ion'
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(8.27, 8.27),
                             num=num, sharex='col', sharey='row')

    # actually plotting the data
    for atupl in ax_tuples:
        
        ax = axes[atupl]
        xind = chan_x[atupl[1]]
        yind = chan_y[atupl[0]]

        energy_x = energy[:, xind]
        energy_y = energy[:, yind]

        ax.plot(
                energy_x[cut_fid], energy_y[cut_fid],
                ls='none', marker='2', zorder=11, color='limegreen',
                label='Fiducial Events'
        )

        ax.plot(
                energy_x[cut_qual], energy_y[cut_qual],
                ls='none', marker='1', zorder=10, color='slateblue',
                label='Quality Events'
        )
            
        ax.plot(
                energy_x, energy_y,
                ls='none', marker=',', zorder=9, color='coral',
                label='All events'
        )

        if xind in run_tree.chan_veto:
            lab = run_tree.chan_label[xind]
            xamp = 2*getattr(noise.sigma0, lab)
            ymin, ymax = energy_y.min(), energy_y.max()
            ax.fill_betweenx([ymin, ymax], -xamp, +xamp, color='lavender')

        if yind in run_tree.chan_veto:
            lab = run_tree.chan_label[yind]
            yamp = 2*getattr(noise.sigma0, lab)
            xmin, xmax = energy_x.min(), energy_x.max()
            ax.fill_between([xmin, xmax], -yamp, +yamp, color='lavender',
                             label='Fiducial selection (2$\sigma$)')
            
        custom_autoscale(ax, energy_x[cut_qual], energy_y[cut_qual])
        
        ax.grid(alpha=0.3)
        
        if atupl == (0,0):
            ax.legend(loc='lower left', framealpha=1,
                      bbox_to_anchor=(1.05, 0.05), borderaxespad=0.,
            )
        
        if atupl[0] == 2:
            ax.set_xlabel(
                    'Energy {} [ADU]'.format(
                            run_tree.chan_label[xind].replace('_', ' ')
                    )
            )
                
        if atupl[1] == 0:
            ax.set_ylabel(
                    'Energy {} [ADU]'.format(
                            run_tree.chan_label[yind].replace('_', ' ')
                    )
            )
    
    fig.text(0.65, 0.98, num,
             horizontalalignment='center',
             verticalalignment='center',
             bbox=dict(facecolor='lime', alpha=0.5))
    
    for tupl in ax_discard:
        fig.delaxes(axes[tupl])
    fig.tight_layout()
    fig.subplots_adjust(hspace=.0, wspace=.0)

#%%
# =============================================================================
# VIRTUAL VS VIRTUAL EV
# =============================================================================

    # recovering data
    energy = trig.energy_ev
    cut_qual = trig.cut.quality
    cut_fid = trig.cut.fiducial
    
    # initializing pseudo-corner plot
    ax_tuples = [(0,0), (1,0), (1,1), (2,0), (2,1), (2,2)]
    ax_discard = [(0, 1), (1, 2), (0, 2)]
    
    
    chan_x = [run_tree.chan_label[ind] for ind in run_tree.chan_collect]
    chan_x.insert(0, 'heat_a')
    chan_y = [run_tree.chan_label[ind] for ind in run_tree.chan_collect]
    chan_y.append('collect')
    
    num = 'VIRTUAL vs VIRTUAL EV'
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(8.27, 8.27),
                             num=num, sharex='col', sharey='row')

    # actually plotting the data
    for atupl in ax_tuples:
        
        ax = axes[atupl]
        xlab = chan_x[atupl[1]]
        ylab = chan_y[atupl[0]]

        energy_x = getattr(energy, xlab)
        energy_y = getattr(energy, ylab)

        ax.plot(
                energy_x[cut_fid], energy_y[cut_fid],
                ls='none', marker='2', zorder=11, color='limegreen',
                label='Fiducial Events'
        )

        ax.plot(
                energy_x[cut_qual], energy_y[cut_qual],
                ls='none', marker='1', zorder=10, color='slateblue',
                label='Quality Events'
        )
            
        ax.plot(
                energy_x, energy_y,
                ls='none', marker=',', zorder=9, color='coral',
                label='All events'
        )
            
        custom_autoscale(ax, energy_x[cut_qual], energy_y[cut_qual])
        
        ax.grid(alpha=0.3)
        
        if atupl == (0,0):
            ax.legend(loc='lower left', framealpha=1,
                      bbox_to_anchor=(1.05, 0.05), borderaxespad=0.,
            )
        
        if atupl[0] == 2:
            ax.set_xlabel(
                    'Energy {} [eV]'.format(
                            xlab.replace('_', ' ')
                    )
            )
                
        if atupl[1] == 0:
            ax.set_ylabel(
                    'Energy {} [eV]'.format(
                            ylab.replace('_', ' ')
                    )
            )
    
    fig.text(0.65, 0.98, num,
             horizontalalignment='center',
             verticalalignment='center',
             bbox=dict(facecolor='lime', alpha=0.5))
    
    for tupl in ax_discard:
        fig.delaxes(axes[tupl])
    fig.tight_layout()
    fig.subplots_adjust(hspace=.0, wspace=.0)
