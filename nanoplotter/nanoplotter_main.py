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


import logging
import sys
import pandas as pd
import numpy as np
from collections import namedtuple
from nanoplotter.plot import Plot
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import seaborn as sns
from pauvre.marginplot import margin_plot
from nanoplotter.timeplots import time_plots
from nanoplotter.spatial_heatmap import spatial_heatmap
from matplotlib import cm
import plotly
import plotly.graph_objs as go


def check_valid_color(color):
    """Check if the color provided by the user is valid.

    If color is invalid the default is returned.
    """
    if color in list(mcolors.CSS4_COLORS.keys()) + ["#4CB391"]:
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
    if colormap in list(cm.cmap_d.keys()):
        logging.info("NanoPlot:  Valid colormap {}.".format(colormap))
        return colormap
    else:
        logging.info("NanoPlot:  Invalid colormap {}, using default.".format(colormap))
        sys.stderr.write("Invalid colormap {}, using default.\n".format(colormap))
        return "Greens"


def check_valid_format(figformat):
    """Check if the specified figure format is valid.

    If format is invalid the default is returned.
    Probably installation-dependent
    """
    fig = plt.figure()
    if figformat in list(fig.canvas.get_supported_filetypes().keys()):
        logging.info("NanoPlot:  valid output format {}".format(figformat))
        return figformat
    else:
        logging.info("NanoPlot:  invalid output format {}".format(figformat))
        sys.stderr.write("Invalid format {}, using default.\n".format(figformat))
        return "png"


def plot_settings(plot_settings, dpi):
    sns.set(**plot_settings)
    mpl.rcParams['savefig.dpi'] = dpi


def scatter(x, y, names, path, plots, color="#4CB391", figformat="png",
            stat=None, log=False, minvalx=0, minvaly=0, title=None,
            plot_settings={}, xmax=None, ymax=None):
    """Create bivariate plots.

    Create four types of bivariate plots of x vs y, containing marginal summaries
    -A scatter plot with histograms on axes
    -A hexagonal binned plot with histograms on axes
    -A kernel density plot with density curves on axes
    -A pauvre-style plot using code from https://github.com/conchoecia/pauvre
    """
    logging.info("NanoPlot:  Creating {} vs {} plots using statistics from {} reads.".format(
        names[0], names[1], x.size))
    if not contains_variance([x, y], names):
        return []
    sns.set(style="ticks", **plot_settings)
    maxvalx = xmax or np.amax(x)
    maxvaly = ymax or np.amax(y)

    plots_made = []

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
            ticks = [10**i for i in range(10) if not 10**i > 10 * (10**maxvalx)]
            plot.ax_joint.set_xticks(np.log10(ticks))
            plot.ax_marg_x.set_xticks(np.log10(ticks))
            plot.ax_joint.set_xticklabels(ticks)
        plt.subplots_adjust(top=0.90)
        plot.fig.suptitle(title or "{} vs {} plot".format(names[0], names[1]), fontsize=25)
        hex_plot.fig = plot
        hex_plot.save(format=figformat)
        plots_made.append(hex_plot)

    sns.set(style="darkgrid", **plot_settings)
    if plots["dot"]:
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
            ticks = [10**i for i in range(10) if not 10**i > 10 * (10**maxvalx)]
            plot.ax_joint.set_xticks(np.log10(ticks))
            plot.ax_marg_x.set_xticks(np.log10(ticks))
            plot.ax_joint.set_xticklabels(ticks)
        plt.subplots_adjust(top=0.90)
        plot.fig.suptitle(title or "{} vs {} plot".format(names[0], names[1]), fontsize=25)
        dot_plot.fig = plot
        dot_plot.save(format=figformat)
        plots_made.append(dot_plot)

    if plots["kde"]:
        idx = np.random.choice(x.index, min(2000, len(x)), replace=False)
        if log:
            kde_plot = Plot(
                path=path + "_loglength_kde." + figformat,
                title="{} vs {} plot using a kernel density estimation "
                      "after log transformation of read lengths".format(names[0], names[1]))
        else:
            kde_plot = Plot(
                path=path + "_kde." + figformat,
                title="{} vs {} plot using a kernel density estimation".format(names[0], names[1]))
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
            ticks = [10**i for i in range(10) if not 10**i > 10 * (10**maxvalx)]
            plot.ax_joint.set_xticks(np.log10(ticks))
            plot.ax_marg_x.set_xticks(np.log10(ticks))
            plot.ax_joint.set_xticklabels(ticks)
        plt.subplots_adjust(top=0.90)
        plot.fig.suptitle(title or "{} vs {} plot".format(names[0], names[1]), fontsize=25)
        kde_plot.fig = plot
        kde_plot.save(format=figformat)
        plots_made.append(kde_plot)

    if plots["pauvre"] and names == ['Read lengths', 'Average read quality'] and log is False:
        pauvre_plot = Plot(
            path=path + "_pauvre." + figformat,
            title="{} vs {} plot using pauvre-style @conchoecia".format(names[0], names[1]))
        sns.set(style="white", **plot_settings)
        margin_plot(df=pd.DataFrame({"length": x, "meanQual": y}),
                    Y_AXES=False,
                    title=title or "Length vs Quality in Pauvre-style",
                    plot_maxlen=None,
                    plot_minlen=0,
                    plot_maxqual=None,
                    plot_minqual=0,
                    lengthbin=None,
                    qualbin=None,
                    BASENAME="whatever",
                    path=pauvre_plot.path,
                    fileform=[figformat],
                    dpi=600,
                    TRANSPARENT=True,
                    QUIET=True)
        plots_made.append(pauvre_plot)
    plt.close("all")
    return plots_made


def contains_variance(arrays, names):
    """
    Make sure both arrays for bivariate ("scatter") plot have a stddev > 0
    """
    for ar, name in zip(arrays, names):
        if np.std(ar) == 0:
            sys.stderr.write(
                "No variation in '{}', skipping bivariate plots.\n".format(name.lower()))
            logging.info("NanoPlot:  No variation in {}, skipping bivariate plot".format(name))
            return False
    else:
        return True


def length_plots(array, name, path, title=None, n50=None, color="#4CB391", figformat="png"):
    """Create histogram of normal and log transformed read lengths."""
    logging.info("NanoPlot:  Creating length plots for {}.".format(name))
    maxvalx = np.amax(array)
    if n50:
        logging.info("NanoPlot:  Using {} reads with read length N50 of {}bp and maximum of {}bp."
                     .format(array.size, n50, maxvalx))
    else:
        logging.info("NanoPlot:  Using {} reads maximum of {}bp.".format(array.size, maxvalx))

    plots = []
    HistType = namedtuple('HistType', 'weight name ylabel')
    for h_type in [HistType(None, "", "Number of reads"),
                   HistType(array, "Weighted ", "Number of bases")]:
        histogram = Plot(
            path=path + h_type.name.replace(" ", "_") + "Histogram" +
            name.replace(' ', '') + "." + figformat,
            title=h_type.name + "Histogram of read lengths")
        ax = sns.distplot(
            a=array,
            kde=False,
            hist=True,
            bins=max(round(int(maxvalx) / 500), 10),
            color=color,
            hist_kws=dict(weights=h_type.weight,
                          edgecolor=color,
                          linewidth=0.2,
                          alpha=0.8))
        if n50:
            plt.axvline(n50)
            plt.annotate('N50', xy=(n50, np.amax([h.get_height() for h in ax.patches])), size=8)
        ax.set(
            xlabel='Read length',
            ylabel=h_type.ylabel,
            title=title or histogram.title)
        plt.ticklabel_format(style='plain', axis='y')
        histogram.fig = ax.get_figure()
        histogram.save(format=figformat)
        plt.close("all")

        log_histogram = Plot(
            path=path + h_type.name.replace(" ", "_") + "LogTransformed_Histogram" +
            name.replace(' ', '') + "." + figformat,
            title=h_type.name + "Histogram of read lengths after log transformation")
        ax = sns.distplot(
            a=np.log10(array),
            kde=False,
            hist=True,
            color=color,
            hist_kws=dict(weights=h_type.weight,
                          edgecolor=color,
                          linewidth=0.2,
                          alpha=0.8))
        ticks = [10**i for i in range(10) if not 10**i > 10 * maxvalx]
        ax.set(
            xticks=np.log10(ticks),
            xticklabels=ticks,
            xlabel='Read length',
            ylabel=h_type.ylabel,
            title=title or log_histogram.title)
        if n50:
            plt.axvline(np.log10(n50))
            plt.annotate('N50', xy=(np.log10(n50), np.amax(
                [h.get_height() for h in ax.patches])), size=8)
        plt.ticklabel_format(style='plain', axis='y')
        log_histogram.fig = ax.get_figure()
        log_histogram.save(format=figformat)
        plt.close("all")
        plots.extend([histogram, log_histogram])
    plots.append(dynamic_histogram(array=array, name=name, path=path, title=title, color=color))
    plots.append(yield_by_minimal_length_plot(array=array,
                                              name=name,
                                              path=path,
                                              title=title,
                                              color=color,
                                              figformat=figformat))
    return plots


def dynamic_histogram(array, name, path, title=None, color="#4CB391"):
    """
    Use plotly to a histogram
    Return html code, but also save as png
    """
    dynhist = Plot(path=path + "Dynamic_Histogram_{}.html".format(name.replace(' ', '_')),
                   title=title or "Dynamic histogram of {}".format(name))
    ylabel = "Number of reads" if len(array) <= 10000 else "Downsampled number of reads"
    dynhist.html, dynhist.fig = plotly_histogram(array=array.sample(min(len(array), 10000)),
                                                 color=color,
                                                 title=dynhist.title,
                                                 xlabel=name,
                                                 ylabel=ylabel)
    dynhist.save()
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
                             title=title)})
    return html, fig


def yield_by_minimal_length_plot(array, name, path,
                                 title=None, color="#4CB391", figformat="png"):
    df = pd.DataFrame(data={"lengths": np.sort(array)[::-1]})
    df["cumyield_gb"] = df["lengths"].cumsum() / 10**9
    yield_by_length = Plot(
        path=path + "Yield_By_Length." + figformat,
        title="Yield by length")
    ax = sns.regplot(
        x='lengths',
        y="cumyield_gb",
        data=df,
        x_ci=None,
        fit_reg=False,
        color=color,
        scatter_kws={"s": 3})
    ax.set(
        xlabel='Read length',
        ylabel='Cumulative yield for minimal length',
        title=title or yield_by_length.title)
    yield_by_length.fig = ax.get_figure()
    yield_by_length.save(format=figformat)
    plt.close("all")
    return yield_by_length


def run_tests():
    import pickle
    df = pickle.load(open("nanotest/sequencing_summary.pickle", "rb"))
    scatter(
        x=df["lengths"],
        y=df["quals"],
        names=['Read lengths', 'Average read quality'],
        path="LengthvsQualityScatterPlot",
        plots={'dot': 1, 'kde': 1, 'hex': 1, 'pauvre': 1},
        plot_settings=dict(font_scale=1))
    time_plots(
        df=df,
        path=".",
        color="#4CB391",
        plot_settings=dict(font_scale=1))
    length_plots(
        array=df["lengths"],
        name="lengths",
        path=".")
    spatial_heatmap(
        array=df["channelIDs"],
        title="Number of reads generated per channel",
        path="ActivityMap_ReadsPerChannel")


if __name__ == "__main__":
    run_tests()
