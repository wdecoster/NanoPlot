# wdecoster
"""
This module provides functions for plotting data extracted from Oxford Nanopore sequencing
reads and alignments, but some of it's functions can also be used for other applications.


FUNCTIONS
* Check if a specified color is a valid matplotlib color
check_valid_color(color)
* Check if a specified output format is valid
checkvalidFormat(format)
* Create a bivariate plot with dots, hexbins and/or kernel density estimates.
Also arguments for specifying axis names, color and xlim/ylim
scatter(x, y, names, path, color, format, plots, stat=None, log=False, minvalx=0, minvaly=0)
* Create cumulative yield plot and evaluate read length and quality over time
timePlots(df, path, color, format)
* Create length distribution histogram and density curve
lengthPlots(array, name, path, n50, color, format, log=False)
* Create flowcell physical layout in numpy array
makeLayout()
* Present the activity (number of reads) per channel on the flowcell as a heatmap
spatialHeatmap(array, title, path, color, format)

"""

import plotly.graph_objs as go
import plotly
import logging
import sys
import os
import pandas as pd
import numpy as np
from nanoplotter.plot import Plot
import plotly.express as px
import plotly.figure_factory as ff
from nanoplotter.spatial_heatmap import spatial_heatmap
from nanoplotter.timeplots import time_plots
import re
import nanoplot.utils as utils


def check_valid_color(color):
    """Check if the color provided by the user is valid.

    If color is invalid the default is returned.
    """
    colors, _ = colors_and_colormaps()
    if color in colors:
        logging.info("NanoPlot:  Valid color {}.".format(color))
        return colors.get(color)

    elif re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', color):
        logging.info("NanoPlot:  Valid color {}.".format(color))
        return color

    else:
        logging.info("NanoPlot:  Invalid color {}, using default.".format(color))
        sys.stderr.write("Invalid color {}, using default.\n".format(color))
        return "#4CB391"


def check_valid_colormap(colormap):
    """Check if the colormap provided by the user is valid.

    If colormap is invalid the default is returned.
    """
    _, colormaps = colors_and_colormaps()
    if colormap in colormaps:
        logging.info("NanoPlot:  Valid colormap {}.".format(colormap))
        return colormap
    else:
        logging.info("NanoPlot:  Invalid colormap {}, using default.".format(colormap))
        sys.stderr.write("Invalid colormap {}, using default.\n".format(colormap))
        return "Greens"


def scatter(x, y, legacy, names, path, plots, color, colormap, figformat, stat=None, log=False, minvalx=0, minvaly=0, title=None, xmax=None, ymax=None):
    """->
    create marginalised scatterplots and KDE plot with marginalized histograms
    -> update from scatter_legacy function to utilise plotly package
    - scatterplot with histogram on both axes
    - kernel density plot with histograms on both axes
    - hexbin not implemented yet
    - pauvre plot temporarily not available
    """
    logging.info("NanoPlot:  Creating {} vs {} plots using statistics from {} reads.".format(
        names[0], names[1], x.size))
    if not contains_variance([x, y], names):
        return []

    plots_made = []
    idx = np.random.choice(x.index, min(10000, len(x)), replace=False)
    maxvalx = xmax or np.amax(x[idx])
    maxvaly = ymax or np.amax(y[idx])

    if plots["dot"]:
        if log:
            dot_plot = Plot(
                path=path + "_loglength_dot.html",
                title=f"{names[0]} vs {names[1]} plot using dots "
                      "after log transformation of read lengths")
        else:
            dot_plot = Plot(
                path=path + "_dot.html",
                title=f"{names[0]} vs {names[1]} plot using dots")

        fig = px.scatter(x=x[idx], y=y[idx], marginal_x="histogram", marginal_y="histogram",
                         range_x=[minvalx, maxvalx], range_y=[minvaly, maxvaly])
        fig.update_traces(marker=dict(color=color))
        fig.update_yaxes(rangemode="tozero")
        fig.update_xaxes(rangemode="tozero")

        fig.update_layout(xaxis_title=names[0],
                          yaxis_title=names[1],
                          title=title or dot_plot.title,
                          title_x=0.5)

        if log:
            ticks = [10 ** i for i in range(10) if not 10 ** i > 10 * (10 ** maxvalx)]
            fig.update_layout(
                xaxis=dict(
                    tickmode='array',
                    tickvals=np.log10(ticks),
                    ticktext=ticks,
                    tickangle=45
                )
            )

        dot_plot.fig = fig
        dot_plot.html = dot_plot.fig.to_html(full_html=False, include_plotlyjs='cdn')
        dot_plot.save(figformat)
        plots_made.append(dot_plot)

    if plots["kde"]:
        if log:
            kde_plot = Plot(
                path=path + "_loglength_kde.html",
                title="{} vs {} plot using a kernel density estimation "
                      "after log transformation of read lengths".format(names[0], names[1]))
        else:
            kde_plot = Plot(
                path=path + "_kde.html",
                title="{} vs {} plot using a kernel density estimation".format(names[0], names[1]))

        col = hex_to_rgb_scale_0_1(color)
        fig = ff.create_2d_density(x[idx], y[idx], point_size=3,
                                   hist_color=col,
                                   point_color=col,
                                   colorscale=colormap, width=1870)

        fig.update_layout(xaxis_title=names[0],
                          yaxis_title=names[1],
                          title=title or kde_plot.title,
                          title_x=0.5,
                          xaxis=dict(tickangle=45))

        if log:
            ticks = [10 ** i for i in range(10) if not 10 ** i > 10 * (10 ** maxvalx)]
            fig.update_layout(
                xaxis=dict(
                    tickmode='array',
                    tickvals=np.log10(ticks),
                    ticktext=ticks,
                    tickangle=45
                )
            )

        kde_plot.fig = fig
        kde_plot.html = kde_plot.fig.to_html(full_html=False, include_plotlyjs='cdn')
        kde_plot.save(figformat)
        plots_made.append(kde_plot)

    if 1 in legacy.values():
        settings, args = utils.get_args()
        plots_made += scatter_legacy(x=x[idx],
                                     y=y[idx],
                                     names=names,
                                     path=path,
                                     plots=legacy,
                                     color=color,
                                     figformat=figformat,
                                     stat=stat,
                                     log=log,
                                     minvalx=minvalx,
                                     minvaly=minvaly,
                                     title=title)
    return plots_made


def scatter_legacy(x, y, names, path, plots, color, figformat,
                   stat=None, log=False, minvalx=0, minvaly=0, title=None,
                   xmax=None, ymax=None):
    """Create bivariate plots.

    Create four types of bivariate plots of x vs y, containing marginal summaries
    -A scatter plot with histograms on axes
    -A hexagonal binned plot with histograms on axes
    -A kernel density plot with density curves on axes
    -A pauvre-style plot using code from https://github.com/conchoecia/pauvre
    """
    try:
        import matplotlib as mpl
        mpl.use('Agg')
        import seaborn as sns
        import matplotlib.pyplot as plt
    except ImportError:
        sys.stderr("need additional modules when running with --legacy")
        return []

    if figformat in ["webp", "json"]:
        figformat = "png"

    logging.info("NanoPlot:  Creating {} vs {} plots using statistics from {} reads (legacy mode).".format(
        names[0], names[1], x.size))
    if not contains_variance([x, y], names):
        return []
    sns.set(style="ticks")
    maxvalx = xmax or np.amax(x)
    maxvaly = ymax or np.amax(y)

    plots_made = []
    path = path + "_legacy"

    if plots["hex"]:
        if log:
            hex_plot = Plot(
                path=path + "_loglength_hex." + figformat,
                title="{} vs {} plot using hexagonal bins "
                      "after log transformation of read lengths".format(names[0], names[1]))
        else:
            hex_plot = Plot(
                path=path + "_hex." + figformat,
                title="{} vs {} plot using hexagonal bins".format(names[0], names[1]))
        plot = sns.jointplot(
            x=x,
            y=y,
            kind="hex",
            color=color,
            stat_func=stat,
            space=0,
            xlim=(minvalx, maxvalx),
            ylim=(minvaly, maxvaly),
            height=10)
        plot.set_axis_labels(names[0], names[1])
        if log:
            ticks = [10 ** i for i in range(10) if not 10 ** i > 10 * (10 ** maxvalx)]
            plot.ax_joint.set_xticks(np.log10(ticks))
            plot.ax_marg_x.set_xticks(np.log10(ticks))
            plot.ax_joint.set_xticklabels(ticks)
        plt.subplots_adjust(top=0.90)
        plot.fig.suptitle(title or "{} vs {} plot".format(names[0], names[1]), fontsize=25)
        hex_plot.fig = plot
        hex_plot.save(figformat)
        plots_made.append(hex_plot)

    sns.set(style="darkgrid")
    if plots["dot"]:
        print("we here")
        if log:
            dot_plot = Plot(
                path=path + "_loglength_dot." + figformat,
                title="{} vs {} plot using dots "
                      "after log transformation of read lengths".format(names[0], names[1]))
        else:
            dot_plot = Plot(
                path=path + "_dot." + figformat,
                title="{} vs {} plot using dots".format(names[0], names[1]))
        plot = sns.jointplot(
            x=x,
            y=y,
            kind="scatter",
            color=color,
            stat_func=stat,
            xlim=(minvalx, maxvalx),
            ylim=(minvaly, maxvaly),
            space=0,
            height=10,
            joint_kws={"s": 1})
        plot.set_axis_labels(names[0], names[1])
        if log:
            ticks = [10 ** i for i in range(10) if not 10 ** i > 10 * (10 ** maxvalx)]
            plot.ax_joint.set_xticks(np.log10(ticks))
            plot.ax_marg_x.set_xticks(np.log10(ticks))
            plot.ax_joint.set_xticklabels(ticks)
        plt.subplots_adjust(top=0.90)
        plot.fig.suptitle(title or "{} vs {} plot".format(names[0], names[1]), fontsize=25)
        dot_plot.fig = plot
        dot_plot.save(figformat)
        plots_made.append(dot_plot)

    if plots["kde"]:
        if len(x) > 2:
            idx = np.random.choice(x.index, min(2000, len(x)), replace=False)
            if log:
                kde_plot = Plot(
                    path=path + "_loglength_kde." + figformat,
                    title="{} vs {} plot using a kernel density estimation "
                          "after log transformation of read lengths".format(names[0], names[1]))
            else:
                kde_plot = Plot(
                    path=path + "_kde." + figformat,
                    title=f"{names[0]} vs {names[1]} plot using a kernel density estimation")
            plot = sns.jointplot(
                x=x[idx],
                y=y[idx],
                kind="kde",
                clip=((0, np.Inf), (0, np.Inf)),
                xlim=(minvalx, maxvalx),
                ylim=(minvaly, maxvaly),
                space=0,
                color=color,
                stat_func=stat,
                shade_lowest=False,
                height=10)
            plot.set_axis_labels(names[0], names[1])
            if log:
                ticks = [10 ** i for i in range(10) if not 10 ** i > 10 * (10 ** maxvalx)]
                plot.ax_joint.set_xticks(np.log10(ticks))
                plot.ax_marg_x.set_xticks(np.log10(ticks))
                plot.ax_joint.set_xticklabels(ticks)
            plt.subplots_adjust(top=0.90)
            plot.fig.suptitle(title or "{} vs {} plot".format(names[0], names[1]), fontsize=25)
            kde_plot.fig = plot
            kde_plot.save(figformat)
            plots_made.append(kde_plot)
        else:
            sys.stderr.write("Not enough observations (reads) to create a kde plot.\n")
            logging.info("NanoPlot: Not enough observations (reads) to create a kde plot")
    plt.close("all")
    return plots_made


# def pauvre_plot():
#     from pauvre.marginplot import margin_plot
#     if plots["pauvre"] and names == ['Read lengths', 'Average read quality'] and log is False:
#         pauvre_plot = Plot(
#             path=path + "_pauvre." + figformat,
#             title="{} vs {} plot using pauvre-style @conchoecia".format(names[0], names[1]))
#         sns.set(style="white")
#         margin_plot(df=pd.DataFrame({"length": x, "meanQual": y}),
#                     Y_AXES=False,
#                     title=title or "Length vs Quality in Pauvre-style",
#                     plot_maxlen=None,
#                     plot_minlen=0,
#                     plot_maxqual=None,
#                     plot_minqual=0,
#                     lengthbin=None,
#                     qualbin=None,
#                     BASENAME="whatever",
#                     path=pauvre_plot.path,
#                     fileform=[figformat],
#                     dpi=600,
#                     TRANSPARENT=True,
#                     QUIET=True)
#         plots_made.append(pauvre_plot)


def contains_variance(arrays, names):
    """
    Make sure both arrays for bivariate ("scatter") plot have a stddev > 0
    """
    for ar, name in zip(arrays, names):
        if np.std(ar) == 0:
            sys.stderr.write(f"No variation in '{name.lower()}', skipping bivariate plots.\n")
            logging.info(f"NanoPlot: No variation in {name}, skipping bivariate plot")
            return False
    else:
        return True


def length_plots(array, name, path, figformat, title=None, n50=None, color="#4CB391"):
    """Create histogram of normal and log transformed read lengths."""
    logging.info("NanoPlot:  Creating length plots for {}.".format(name))
    maxvalx = np.amax(array)
    if n50:
        logging.info("NanoPlot: Using {} reads with read length N50 of {}bp and maximum of {}bp."
                     .format(array.size, n50, maxvalx))
    else:
        logging.info(f"NanoPlot:  Using {array.size} reads maximum of {maxvalx}bp.")

    plots = []

    HistType = [{'weight': array, 'name': 'Weighted', 'ylabel': 'Number of reads'},
                {'weight': None, 'name': 'Non weighted', 'ylabel': 'Number of reads'}]

    for h_type in HistType:
        histogram = Plot(
            path=path + h_type["name"].replace(" ", "_") + "Histogram" +
                 name.replace(' ', '') + ".html",
            title=f"{h_type['name']} histogram of read lengths")

        hist, bin_edges = np.histogram(array,
                                       bins=max(round(int(maxvalx) / 500), 10),
                                       weights=h_type["weight"])

        fig = go.Figure()

        fig.add_trace(go.Bar(x=bin_edges[1:],
                             y=hist,
                             marker_color=color))

        if n50:
            fig.add_vline(n50)
            fig.add_annotation(text='N50', x=n50, y=0.95)
            fig.update_annotations(font_size=8)

        fig.update_layout(xaxis_title='Read length',
                          yaxis_title=h_type["ylabel"],
                          title=title or histogram.title,
                          title_x=0.5)

        histogram.fig = fig
        histogram.html = histogram.fig.to_html(full_html=False, include_plotlyjs='cdn')
        histogram.save(figformat)

        log_histogram = Plot(
            path=path + h_type["name"].replace(" ", "_") + "LogTransformed_Histogram" +
                 name.replace(' ', '') + ".html",
            title=h_type["name"] + " histogram of read lengths after log transformation")

        if h_type["weight"] is None:
            hist_log, bin_edges_log = np.histogram(np.log10(array),
                                                   bins=max(round(int(maxvalx) / 500), 10),
                                                   weights=h_type["weight"])

        else:
            hist_log, bin_edges_log = np.histogram(np.log10(array),
                                                   bins=max(round(int(maxvalx) / 500), 10),
                                                   weights=np.log10(h_type["weight"]))

        fig = go.Figure()
        fig.add_trace(go.Bar(x=bin_edges_log[1:],
                             y=hist_log,
                             marker_color=color))

        ticks = [10 ** i for i in range(10) if not 10 ** i > 10 * maxvalx]

        fig.update_layout(
            xaxis=dict(
                tickmode='array',
                tickvals=np.log10(ticks),
                ticktext=ticks),
            xaxis_title='Read length',
            yaxis_title=h_type["ylabel"],
            title=title or log_histogram.title,
            title_x=0.5)

        if n50:
            fig.add_vline(np.log10(n50))
            fig.add_annotation(text='N50', x=np.log10(n50), y=0.95)
            fig.update_annotations(font_size=8)

        log_histogram.fig = fig
        log_histogram.html = log_histogram.fig.to_html(full_html=False, include_plotlyjs='cdn')
        log_histogram.save(figformat)

        plots.extend([histogram, log_histogram])

    plots.append(yield_by_minimal_length_plot(array=array,
                                                name=name,
                                                path=path,
                                                title=title,
                                                color=color,
                                                figformat=figformat))

    return plots


def dynamic_histogram(array, name, path, figformat, title=None, color="#4CB391"):
    """
    Use plotly to a histogram
    Return html code, but also save as png
    """
    dynhist = Plot(
        path=path + f"Dynamic_Histogram_{name[0].lower() + name[1:].replace(' ', '_')}.html",
        title="Dynamic histogram of {}".format(name[0].lower() + name[1:]))
    ylabel = "Number of reads" if len(array) <= 10000 else "Downsampled number of reads"
    dynhist.html, dynhist.fig = plotly_histogram(array=array.sample(min(len(array), 10000)),
                                                 color=color,
                                                 title=title or dynhist.title,
                                                 xlabel=name,
                                                 ylabel=ylabel)
    dynhist.save(figformat)
    return dynhist


def plotly_histogram(array, color="#4CB391", title=None, xlabel=None, ylabel=None):
    data = [go.Histogram(x=array,
                         opacity=0.4,
                         marker=dict(color=color))]
    html = plotly.offline.plot(
        {"data": data,
         "layout": go.Layout(barmode='overlay',
                             title=title,
                             yaxis_title=ylabel,
                             xaxis_title=xlabel)},
        output_type="div",
        show_link=False)
    fig = go.Figure(
        {"data": data,
         "layout": go.Layout(barmode='overlay',
                             title=title,
                             title_x=0.5)})
    return html, fig


def yield_by_minimal_length_plot(array, name, path, figformat, title=None, color="#4CB391"):
    df = pd.DataFrame(data={"lengths": np.sort(array)[::-1]})
    df["cumyield_gb"] = df["lengths"].cumsum() / 10 ** 9
    idx = np.random.choice(array.index, min(10000, len(array)), replace=False)

    yield_by_length = Plot(
        path=path + "Yield_By_Length.html",
        title="Yield by length")

    fig = px.scatter(df,x=df.reindex(idx)["lengths"], y=df.reindex(idx)["cumyield_gb"])
    fig.update_traces(marker=dict(color=color))
    fig.update_layout(xaxis_title='Read length',
                      yaxis_title='Cumulative yield for minimal length [Gb]',
                      title=title or yield_by_length.title,
                      title_x=0.5)

    yield_by_length.fig = fig
    yield_by_length.html = yield_by_length.fig.to_html(full_html=False, include_plotlyjs='cdn')
    yield_by_length.save(figformat)

    return yield_by_length


def colors_and_colormaps():
    colormaps = ('Greys,YlGnBu,Greens,YlOrRd,Bluered,RdBu,Reds,Blues,Picnic,Rainbow,Portland,Jet,'
                 'Hot,Blackbody,Earth,Electric,Viridis,Cividis').split(',')
    parent_directory = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
    colours = open(os.path.join(parent_directory, "extra/color_options_hex.txt"))
    col_hex = {}

    for line in colours:
        key, value = line.split(",")
        col_hex[key] = value.strip()

    return col_hex, colormaps


def hex_to_rgb_scale_0_1(hexcolor):
    color = hexcolor.lstrip("#")
    RGB_color = tuple(int(color[x:x + 2], 16) for x in (0, 2, 4))

    RGB_color = [x / 255 for x in RGB_color]

    return tuple(RGB_color)


def run_tests():
    import pickle
    df = pickle.load(open("nanotest/sequencing_summary.pickle", "rb"))
    scatter(
        x=df["lengths"],
        y=df["quals"],
        names=['Read lengths', 'Average read quality'],
        path="LengthvsQualityScatterPlot",
        plots={'dot': 1, 'kde': 1})
    time_plots(
        df=df,
        path="./",
        color="#4CB391")
    length_plots(
        array=df["lengths"],
        name="lengths",
        path="./")
    spatial_heatmap(
        array=df["channelIDs"],
        title="Number of reads generated per channel",
        path="ActivityMap_ReadsPerChannel")


if __name__ == "__main__":
    run_tests()
