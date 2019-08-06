#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
@author: misiak

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe


from core_classes import Analysis, Artifact
from plot_addon import LegendTitle, custom_autoscale, ax_hist
from model_spectrum import double_norm, fid_mixture
from stats_addon import cdf_calc, custom_bin_edges

from spec_classes import Analysis_red80
#%%

# first command
plt.close('all')
plt.rcParams['text.usetex']=True
style = [pe.Normal(), pe.withStroke(foreground='k', linewidth=3)]

run = 'tg25l019'
ana = Analysis_red80(run)

trig = ana.all.trig
noise = ana.all.noise
run_tree = ana.all.run_tree

etypes = [trig, noise]

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
    
    # handler_map + LegendTitle allow for subtitle in legend
    ax.legend(
            handles, labels, loc=2, framealpha=1,
            bbox_to_anchor=(1.05, 1), borderaxespad=0.,
            handler_map={str: LegendTitle()}
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
                cdf = ana.model.dist.cdf(xrange, *popt)
                normalization = getattr(trig,
                                        'nsamples_{}'.format(
                                                ana.calibration_peak.cut_type
                                        ))
                pdf_norm = pdf * normalization * (bin_edges[1] - bin_edges[0])
                
                ax.autoscale(False)
                ax.plot(xrange, pdf_norm,
                        ls='--', color='yellow',
                        label='fit')
                
                ax.twinx().plot(xrange, cdf,
                        ls='-.', color='yellow',
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
