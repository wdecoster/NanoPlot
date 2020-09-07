import sys
import logging
from nanoplotter.plot import Plot
from datetime import timedelta
import seaborn as sns
import matplotlib.pyplot as plt
from math import ceil
import pandas as pd
import numpy as np


def check_valid_time_and_sort(df, timescol, days=5, warning=True):
    """Check if the data contains reads created within the same `days` timeframe.

    if not, print warning and only return part of the data which is within `days` days
    Resetting the index twice to get also an "index" column for plotting the cum_yield_reads plot
    """
    timediff = (df[timescol].max() - df[timescol].min()).days
    if timediff < days:
        return df.sort_values(timescol).reset_index(drop=True).reset_index()
    else:
        if warning:
            sys.stderr.write(
                "\nWarning: data generated is from more than {} days.\n".format(str(days)))
            sys.stderr.write("Likely this indicates you are combining multiple runs.\n")
            sys.stderr.write(
                "Plots based on time are invalid and therefore truncated to first {} days.\n\n"
                .format(str(days)))
            logging.warning("Time plots truncated to first {} days: invalid timespan: {} days"
                            .format(str(days), str(timediff)))
        return df[df[timescol] < timedelta(days=days)] \
            .sort_values(timescol) \
            .reset_index(drop=True) \
            .reset_index()


def time_plots(df, path, title=None, color="#4CB391", figformat="png",
               log_length=False, plot_settings=None):
    """Making plots of time vs read length, time vs quality and cumulative yield."""
    dfs = check_valid_time_and_sort(df, "start_time")
    logging.info("Nanoplotter: Creating timeplots using {} reads.".format(len(dfs)))
    cumyields = cumulative_yield(dfs=dfs.set_index("start_time"),
                                 path=path,
                                 figformat=figformat,
                                 title=title,
                                 color=color)
    reads_pores_over_time = plot_over_time(dfs=dfs.set_index("start_time"),
                                           path=path,
                                           figformat=figformat,
                                           title=title,
                                           color=color)
    violins = violin_plots_over_time(dfs=dfs,
                                     path=path,
                                     figformat=figformat,
                                     title=title,
                                     log_length=log_length,
                                     plot_settings=plot_settings)
    return cumyields + reads_pores_over_time + violins


def violin_plots_over_time(dfs, path, figformat, title,
                           log_length=False, plot_settings=None):
    dfs['timebin'] = add_time_bins(dfs)
    plots = []
    plots.append(length_over_time(dfs=dfs,
                                  path=path,
                                  figformat=figformat,
                                  title=title,
                                  log_length=log_length,
                                  plot_settings=plot_settings))
    if "quals" in dfs:
        plots.append(quality_over_time(dfs=dfs,
                                       path=path,
                                       figformat=figformat,
                                       title=title,
                                       plot_settings=plot_settings))
    if "duration" in dfs:
        plots.append(sequencing_speed_over_time(dfs=dfs,
                                                path=path,
                                                figformat=figformat,
                                                title=title,
                                                plot_settings=plot_settings))
    return plots


def length_over_time(dfs, path, figformat, title, log_length=False, plot_settings={}):
    if log_length:
        time_length = Plot(path=path + "TimeLogLengthViolinPlot." + figformat,
                           title="Violin plot of log read lengths over time")
    else:
        time_length = Plot(path=path + "TimeLengthViolinPlot." + figformat,
                           title="Violin plot of read lengths over time")
    sns.set(style="white", **plot_settings)
    if log_length:
        length_column = "log_lengths"
    else:
        length_column = "lengths"

    if "length_filter" in dfs:  # produced by NanoPlot filtering of too long reads
        temp_dfs = dfs[dfs["length_filter"]]
    else:
        temp_dfs = dfs

    ax = sns.violinplot(x="timebin",
                        y=length_column,
                        data=temp_dfs,
                        inner=None,
                        cut=0,
                        linewidth=0)
    ax.set(xlabel='Interval (hours)',
           ylabel="Read length",
           title=title or time_length.title)
    if log_length:
        ticks = [10**i for i in range(10) if not 10**i > 10 * np.amax(dfs["lengths"])]
        ax.set(yticks=np.log10(ticks),
               yticklabels=ticks)
    plt.xticks(rotation=45, ha='center', fontsize=8)
    time_length.fig = ax.get_figure()
    time_length.save(format=figformat)
    plt.close("all")
    return time_length


def quality_over_time(dfs, path, figformat, title, plot_settings={}):
    time_qual = Plot(path=path + "TimeQualityViolinPlot." + figformat,
                     title="Violin plot of quality over time")
    sns.set(style="white", **plot_settings)
    ax = sns.violinplot(x="timebin",
                        y="quals",
                        data=dfs,
                        inner=None,
                        cut=0,
                        linewidth=0)
    ax.set(xlabel='Interval (hours)',
           ylabel="Basecall quality",
           title=title or time_qual.title)
    plt.xticks(rotation=45, ha='center', fontsize=8)
    time_qual.fig = ax.get_figure()
    time_qual.save(format=figformat)
    plt.close("all")
    return time_qual


def sequencing_speed_over_time(dfs, path, figformat, title, plot_settings={}):
    time_duration = Plot(path=path + "TimeSequencingSpeed_ViolinPlot." + figformat,
                         title="Violin plot of sequencing speed over time")
    sns.set(style="white", **plot_settings)
    if "timebin" not in dfs:
        dfs['timebin'] = add_time_bins(dfs)
    mask = dfs['duration'] != 0
    ax = sns.violinplot(x=dfs.loc[mask, "timebin"],
                        y=dfs.loc[mask, "lengths"] / dfs.loc[mask, "duration"],
                        inner=None,
                        cut=0,
                        linewidth=0)
    ax.set(xlabel='Interval (hours)',
           ylabel="Sequencing speed (nucleotides/second)",
           title=title or time_duration.title)
    plt.xticks(rotation=45, ha='center', fontsize=8)
    time_duration.fig = ax.get_figure()
    time_duration.save(format=figformat)
    plt.close("all")
    return time_duration


def add_time_bins(dfs, bin_length=3):
    maxtime = dfs["start_time"].max().total_seconds()
    labels = [str(i) + "-" + str(i + bin_length)
              for i in range(0, 168, bin_length) if not i >= (maxtime / 3600)]
    return pd.cut(x=dfs["start_time"],
                  bins=ceil((maxtime / 3600) / bin_length),
                  labels=labels)


def plot_over_time(dfs, path, figformat, title, color):
    num_reads = Plot(path=path + "NumberOfReads_Over_Time." + figformat,
                     title="Number of reads over time")
    s = dfs.loc[:, "lengths"].resample('10T').count()
    ax = sns.regplot(x=s.index.total_seconds() / 3600,
                     y=s,
                     x_ci=None,
                     fit_reg=False,
                     color=color,
                     scatter_kws={"s": 3})
    ax.set(xlabel='Run time (hours)',
           ylabel='Number of reads per 10 minutes',
           title=title or num_reads.title)
    num_reads.fig = ax.get_figure()
    num_reads.save(format=figformat)
    plt.close("all")
    plots = [num_reads]

    if "channelIDs" in dfs:
        pores_over_time = Plot(path=path + "ActivePores_Over_Time." + figformat,
                               title="Number of active pores over time")
        s = dfs.loc[:, "channelIDs"].resample('10T').nunique()
        ax = sns.regplot(x=s.index.total_seconds() / 3600,
                         y=s,
                         x_ci=None,
                         fit_reg=False,
                         color=color,
                         scatter_kws={"s": 3})
        ax.set(xlabel='Run time (hours)',
               ylabel='Active pores per 10 minutes',
               title=title or pores_over_time.title)
        pores_over_time.fig = ax.get_figure()
        pores_over_time.save(format=figformat)
        plt.close("all")
        plots.append(pores_over_time)
    return plots


def cumulative_yield(dfs, path, figformat, title, color):
    cum_yield_gb = Plot(path=path + "CumulativeYieldPlot_Gigabases." + figformat,
                        title="Cumulative yield")
    s = dfs.loc[:, "lengths"].cumsum().resample('1T').max() / 1e9
    ax = sns.regplot(x=s.index.total_seconds() / 3600,
                     y=s,
                     x_ci=None,
                     fit_reg=False,
                     color=color,
                     scatter_kws={"s": 3})
    ax.set(xlabel='Run time (hours)',
           ylabel='Cumulative yield in gigabase',
           title=title or cum_yield_gb.title)
    cum_yield_gb.fig = ax.get_figure()
    cum_yield_gb.save(format=figformat)
    plt.close("all")

    cum_yield_reads = Plot(path=path + "CumulativeYieldPlot_NumberOfReads." + figformat,
                           title="Cumulative yield")
    s = dfs.loc[:, "lengths"].resample('10T').count().cumsum()
    ax = sns.regplot(x=s.index.total_seconds() / 3600,
                     y=s,
                     x_ci=None,
                     fit_reg=False,
                     color=color,
                     scatter_kws={"s": 3})
    ax.set(xlabel='Run time (hours)',
           ylabel='Cumulative yield in number of reads',
           title=title or cum_yield_reads.title)
    cum_yield_reads.fig = ax.get_figure()
    cum_yield_reads.save(format=figformat)
    plt.close("all")
    return [cum_yield_gb, cum_yield_reads]
