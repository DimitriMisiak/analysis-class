#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Some functions and classes for graphic representation.

@author: misiak
"""

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np

from stats_addon import cdf_calc

def lighten_color(color, amount=0.5):
    """
    Credits to: 
    Ian Hincks ihincks
    https://gist.github.com/ihincks/6a420b599f43fcd7dbd79d56798c4e5a
    
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.
    
    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0],1-amount * (1-c[1]),c[2])


# subtitles in legend
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


def ax_hist(axis, bin_edges, data_array, lab, color='slateblue'):
    """ Draw pretty histogramm and cdf in given axis.
    Return bin_array, hist_array, cdf_array.    
    """
    c_dark = lighten_color(color, 1.5)
    c_light = lighten_color(color, 0.8)
    
    style = [pe.Normal(), pe.withStroke(foreground='k', linewidth=3)]

    bin_array = bin_edges[1:] - (bin_edges[1]-bin_edges[0])/2

    data_hist, _ = np.histogram(data_array, bins=bin_edges)
    
    data_sorted, cdf_array = cdf_calc(data_array)
    
    hist_line, = axis.plot(bin_array, data_hist, ls='steps-mid',
                           color=c_dark)

    axis.fill_between(bin_array, data_hist, label=lab,
                      color=c_light, step='mid')

    a0 = axis.twinx()
    a0.set_ylabel('CDF', color='grey')
    a0.tick_params(axis='y', labelcolor='grey')
    
    cdf_line, = a0.plot(data_sorted, cdf_array,
                        ls='steps', color=color, path_effects=style)
    
    axis.grid(True)
    axis.set_ylabel('Counts Events {}'.format(lab), color='k')
    axis.tick_params(axis='y', labelcolor='k')
    axis.set_xlabel('Energy [ADU]')
    
    axis.legend(loc=2)
    
    axis.set_yscale('log')
    axis.set_xlim(bin_edges[0], bin_edges[-1])
    
    return bin_array, data_hist, cdf_array

