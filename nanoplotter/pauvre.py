#!/usr/bin/env python
# -*- coding: utf-8 -*-

# pauvre - just a pore PhD student's plotting package
# Copyright (c) 2016-2017 Darrin T. Schultz. All rights reserved.
#
# This file is part of pauvre.
#
# pauvre is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pauvre is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pauvre.  If not, see <http://www.gnu.org/licenses/>.

# 20191204 This code was modified by wdecoster to be specifically tailored to NanoPlot


import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rcParams
import numpy as np
import logging

# logging
logger = logging.getLogger('pauvre')


def generate_panel(panel_left, panel_bottom, panel_width, panel_height,
                   axis_tick_param='both', which_tick_param='both',
                   bottom_tick_param='on', label_bottom_tick_param='on',
                   left_tick_param='on', label_left_tick_param='on',
                   right_tick_param='off', label_right_tick_param='off',
                   top_tick_param='off', label_top_tick_param='off'):
    """
        Setting default panel tick parameters. Some of these are the defaults
        for matplotlib anyway, but specifying them for readability. Here are
        options and defaults for the parameters used below:

        axis : {'x', 'y', 'both'}; which axis to modify; default = 'both'
        which : {'major', 'minor', 'both'}; which ticks to modify;
                default = 'major'
        bottom, top, left, right : bool or {'on', 'off'}; ticks on or off;
        labelbottom, labeltop, labelleft, labelright : bool or {'on', 'off'}
     """

    # create the panel
    panel_rectangle = [panel_left, panel_bottom, panel_width, panel_height]
    panel = plt.axes(panel_rectangle)

    # Set tick parameters
    panel.tick_params(axis=axis_tick_param,
                      which=which_tick_param,
                      bottom=bottom_tick_param,
                      labelbottom=label_bottom_tick_param,
                      left=left_tick_param,
                      labelleft=label_left_tick_param,
                      right=right_tick_param,
                      labelright=label_right_tick_param,
                      top=top_tick_param,
                      labeltop=label_top_tick_param)

    return panel


def _generate_histogram_bin_patches(panel, bins, bin_values, horizontal=True):
    """This helper method generates the histogram that is added to the panel.

    In this case, horizontal = True applies to the mean quality histogram.
    So, horizontal = False only applies to the length histogram.
    """
    if horizontal:
        for step in np.arange(0, len(bin_values), 1):
            left = bins[step]
            bottom = 0
            width = bins[step + 1] - bins[step]
            height = bin_values[step]
            hist_rectangle = mplpatches.Rectangle((left, bottom), width, height,
                                                  linewidth=0.0,
                                                  facecolor=(0.5, 0.5, 0.5),
                                                  edgecolor=(0, 0, 0))
            panel.add_patch(hist_rectangle)
    else:
        for step in np.arange(0, len(bin_values), 1):
            left = 0
            bottom = bins[step]
            width = bin_values[step]
            height = bins[step + 1] - bins[step]

            hist_rectangle = mplpatches.Rectangle((left, bottom), width, height,
                                                  linewidth=0.0,
                                                  facecolor=(0.5, 0.5, 0.5),
                                                  edgecolor=(0, 0, 0))
            panel.add_patch(hist_rectangle)


def generate_histogram(panel, data_list, max_plot_length, min_plot_length,
                       bin_interval, hist_horizontal=True,
                       left_spine=True, bottom_spine=True,
                       top_spine=False, right_spine=False, x_label=None,
                       y_label=None):

    bins = np.arange(0, max_plot_length, bin_interval or 1)

    bin_values, bins2 = np.histogram(data_list, bins)

    # hist_horizontal is used for quality
    if hist_horizontal:
        panel.set_xlim([min_plot_length, max_plot_length])
        panel.set_ylim([0, max(bin_values * 1.1)])
    # and hist_horizontal == Fale is for read length
    else:
        panel.set_xlim([0, max(bin_values * 1.1)])
        panel.set_ylim([min_plot_length, max_plot_length])

    # Generate histogram bin patches, depending on whether we're plotting
    # vertically or horizontally
    _generate_histogram_bin_patches(panel, bins, bin_values, hist_horizontal)

    panel.spines['left'].set_visible(left_spine)
    panel.spines['bottom'].set_visible(bottom_spine)
    panel.spines['top'].set_visible(top_spine)
    panel.spines['right'].set_visible(right_spine)

    if y_label is not None:
        panel.set_ylabel(y_label)
    if x_label is not None:
        panel.set_xlabel(x_label)


def generate_heat_map(panel, data_frame, min_plot_length, min_plot_qual,
                      max_plot_length, max_plot_qual, color, **kwargs):
    panel.set_xlim([min_plot_qual, max_plot_qual])
    panel.set_ylim([min_plot_length, max_plot_length])

    hex_this = data_frame.query('length<{} and meanQual<{}'.format(
        max_plot_length, max_plot_qual))

    # This single line controls plotting the hex bins in the panel
    hex_vals = panel.hexbin(hex_this['meanQual'], hex_this['length'], gridsize=49,
                            linewidths=0.0, cmap=color)
    for each in panel.spines:
        panel.spines[each].set_visible(False)

    counts = hex_vals.get_array()
    return counts


def generate_legend(panel, counts, color):

    # completely custom for more control
    panel.set_xlim([0, 1])
    panel.set_ylim([0, 1000])
    panel.set_yticks([int(x) for x in np.linspace(0, 1000, 6)])
    panel.set_yticklabels([int(x) for x in np.linspace(0, max(counts), 6)])
    for i in np.arange(0, 1001, 1):
        rgba = color(i / 1001)
        alpha = rgba[-1]
        facec = rgba[0:3]
        hist_rectangle = mplpatches.Rectangle((0, i), 1, 1,
                                              linewidth=0.0,
                                              facecolor=facec,
                                              edgecolor=(0, 0, 0),
                                              alpha=alpha)
        panel.add_patch(hist_rectangle)
    panel.spines['top'].set_visible(False)
    panel.spines['left'].set_visible(False)
    panel.spines['bottom'].set_visible(False)
    panel.yaxis.set_label_position("right")
    panel.set_ylabel('Number of Reads')


def update_rcParams():
    # This mpl style is from the UCSC BME163 class.
    rcParams.update({
        'font.size': 8.0,
        'font.sans-serif': 'Arial',
        'xtick.major.size': 2,
        'xtick.major.width': 0.75,
        'xtick.labelsize': 8.0,
        'xtick.direction': 'out',
        'ytick.major.size': 2,
        'ytick.major.width': 0.75,
        'ytick.labelsize': 8.0,
        'ytick.direction': 'out',
        'xtick.major.pad': 2,
        'xtick.minor.pad': 2,
        'ytick.major.pad': 2,
        'ytick.minor.pad': 2,
        'savefig.dpi': 601,
        'axes.linewidth': 0.75,
        'text.usetex': False})


def margin_plot(df, **kwargs):
    update_rcParams()

    # 250, 231, 34 light yellow
    # 67, 1, 85
    # R=np.linspace(65/255,1,101)
    # G=np.linspace(0/255, 231/255, 101)
    # B=np.linspace(85/255, 34/255, 101)
    # R=65/255, G=0/255, B=85/255
    Rf = 65 / 255
    Bf = 85 / 255
    pdict = {'red': ((0.0, Rf, Rf),
                     (1.0, Rf, Rf)),
             'green': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
             'blue': ((0.0, Bf, Bf),
                      (1.0, Bf, Bf)),
             'alpha': ((0.0, 0.0, 0.0),
                       (1.0, 1.0, 1.0))
             }
    # Now we will use this example to illustrate 3 ways of
    # handling custom colormaps.
    # First, the most direct and explicit:
    purple1 = LinearSegmentedColormap('Purple1', pdict)

    # set the figure dimensions
    fig_width = 1.61 * 3
    fig_height = 1 * 3
    fig = plt.figure(figsize=(fig_width, fig_height))

    # set the panel dimensions
    heat_map_panel_width = fig_width * 0.5
    heat_map_panel_height = heat_map_panel_width * 0.62

    # find the margins to center the panel in figure
    fig_left_margin = fig_bottom_margin = (1 / 6)

    # lengthPanel
    length_panel_width = (1 / 8)

    # the color Bar parameters
    legend_panel_width = (1 / 24)

    # define padding
    h_padding = 0.02
    v_padding = 0.05

    # Set whether to include y-axes in histograms
    length_bottom_spine = False
    length_bottom_tick = 'off'
    length_bottom_label = 'off'
    qual_left_spine = False
    qual_left_tick = 'off'
    qual_left_label = 'off'
    qual_y_label = None

    panels = []

    # Quality histogram panel
    qual_panel_left = fig_left_margin + length_panel_width + h_padding
    qual_panel_width = heat_map_panel_width / fig_width
    qual_panel_height = length_panel_width * fig_width / fig_height
    qual_panel = generate_panel(qual_panel_left,
                                fig_bottom_margin,
                                qual_panel_width,
                                qual_panel_height,
                                left_tick_param=qual_left_tick,
                                label_left_tick_param=qual_left_label)
    panels.append(qual_panel)

    # Length histogram panel
    length_panel_bottom = fig_bottom_margin + qual_panel_height + v_padding
    length_panel_height = heat_map_panel_height / fig_height
    length_panel = generate_panel(fig_left_margin,
                                  length_panel_bottom,
                                  length_panel_width,
                                  length_panel_height,
                                  bottom_tick_param=length_bottom_tick,
                                  label_bottom_tick_param=length_bottom_label)
    panels.append(length_panel)

    # Heat map panel
    heat_map_panel_left = fig_left_margin + length_panel_width + h_padding
    heat_map_panel_bottom = fig_bottom_margin + qual_panel_height + v_padding

    heat_map_panel = generate_panel(heat_map_panel_left,
                                    heat_map_panel_bottom,
                                    heat_map_panel_width / fig_width,
                                    heat_map_panel_height / fig_height,
                                    bottom_tick_param='off',
                                    label_bottom_tick_param='off',
                                    left_tick_param='off',
                                    label_left_tick_param='off')
    panels.append(heat_map_panel)
    heat_map_panel.set_title(kwargs["title"])

    # Legend panel
    legend_panel_left = fig_left_margin + length_panel_width + \
        heat_map_panel_width / fig_width + h_padding
    legend_panel_bottom = fig_bottom_margin + qual_panel_height + v_padding
    legend_panel_height = heat_map_panel_height / fig_height
    legend_panel = generate_panel(legend_panel_left, legend_panel_bottom,
                                  legend_panel_width, legend_panel_height,
                                  bottom_tick_param='off',
                                  label_bottom_tick_param='off',
                                  left_tick_param='off',
                                  label_left_tick_param='off',
                                  right_tick_param='on',
                                  label_right_tick_param='on')
    panels.append(legend_panel)

    # Set min and max viewing window for length
    max_plot_length = int(np.percentile(df['length'], 99))
    min_plot_length = 0

    # Set length bin sizes
    length_bin_interval = int(max_plot_length / 80)

    # length_bins = np.arange(0, max_plot_length, length_bin_interval)

    # Set max and min viewing window for quality
    max_plot_qual = max(np.ceil(df['meanQual']))
    min_plot_qual = 0

    # Set qual bin sizes
    qual_bin_interval = max_plot_qual / 85

    # Generate length histogram
    generate_histogram(length_panel, df['length'], max_plot_length, min_plot_length,
                       length_bin_interval, hist_horizontal=False,
                       y_label='Read Length', bottom_spine=length_bottom_spine)

    # Generate quality histogram
    generate_histogram(qual_panel, df['meanQual'], max_plot_qual, min_plot_qual,
                       qual_bin_interval, x_label='Phred Quality',
                       y_label=qual_y_label, left_spine=qual_left_spine)

    # Generate heat map
    counts = generate_heat_map(heat_map_panel, df, min_plot_length, min_plot_qual,
                               max_plot_length, max_plot_qual, purple1)

    # Generate legend
    generate_legend(legend_panel, counts, purple1)

    return fig
