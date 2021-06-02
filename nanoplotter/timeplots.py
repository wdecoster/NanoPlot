import sys
import logging
from nanoplotter.plot import Plot
from datetime import timedelta
from math import ceil
import pandas as pd
import numpy as np
import plotly.graph_objs as go
import plotly.express as px


def check_valid_time_and_sort(df, timescol="start_time", days=5, warning=True):
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


def time_plots(df, subsampled_df, path, figformat, title=None, color="#4CB391", log_length=False):
    """Making plots of time vs read length, time vs quality and cumulative yield."""

    logging.info(f"Nanoplotter: Creating timeplots using {len(df)} (full) or "
                 f"{len(subsampled_df)} (subsampled dataset) reads.")
    dfs = check_valid_time_and_sort(df)
    cumyields = cumulative_yield(dfs=dfs.set_index("start_time"),
                                 path=path,
                                 title=title,
                                 color=color,
                                 figformat=figformat)
    reads_pores_over_time = plot_over_time(dfs=dfs.set_index("start_time"),
                                           path=path,
                                           title=title,
                                           color=color,
                                           figformat=figformat)
    violins = violin_plots_over_time(dfs=check_valid_time_and_sort(subsampled_df),
                                     path=path,
                                     title=title,
                                     log_length=log_length,
                                     color=color,
                                     figformat=figformat)
    return cumyields + reads_pores_over_time + violins


def violin_plots_over_time(dfs, path, title, figformat, log_length=False, color="#4CB391"):

    dfs['timebin'] = add_time_bins(dfs)
    plots = []

    dfs.sort_values("timebin")

    plots.append(length_over_time(dfs=dfs,
                                  path=path,
                                  title=title,
                                  log_length=log_length,
                                  color=color,
                                  figformat=figformat))
    if "quals" in dfs:
        plots.append(quality_over_time(dfs=dfs,
                                       path=path,
                                       title=title,
                                       color=color,
                                       figformat=figformat))
    if "duration" in dfs:
        plots.append(sequencing_speed_over_time(dfs=dfs,
                                                path=path,
                                                title=title,
                                                color=color,
                                                figformat=figformat))
    return plots


def length_over_time(dfs, path, title, figformat, log_length=False, color="#4CB391"):
    if log_length:
        time_length = Plot(path=path + "TimeLogLengthViolinPlot.html",
                           title="Violin plot of log read lengths over time")
    else:
        time_length = Plot(path=path + "TimeLengthViolinPlot.html",
                           title="Violin plot of read lengths over time")

    length_column = "log_lengths" if log_length else "lengths"

    if "length_filter" in dfs:  # produced by NanoPlot filtering of too long reads
        temp_dfs = dfs[dfs["length_filter"]]
    else:
        temp_dfs = dfs

    fig = go.Figure()

    fig.add_trace(go.Violin(y=temp_dfs[length_column],
                            x=temp_dfs["timebin"],
                            points=False, spanmode="hard",
                            line_color='black', line_width=1.5,
                            fillcolor=color, opacity=0.8))
    fig.update_layout(xaxis_title='Interval (hours)',
                      yaxis_title='Read length',
                      title=title or time_length.title,
                      title_x=0.5)

    if log_length:
        ticks = [10**i for i in range(10) if not 10**i > 10 * np.amax(dfs["lengths"])]
        fig.update_layout(
            yaxis=dict(
                tickmode='array',
                tickvals=np.log10(ticks),
                ticktext=ticks
            )
        )

    fig.update_yaxes(tickangle=45)

    time_length.fig = fig
    time_length.html = time_length.fig.to_html(full_html=False, include_plotlyjs='cdn')
    time_length.save(figformat)

    return time_length


def quality_over_time(dfs, path, figformat, title=None, color="#4CB391"):
    time_qual = Plot(path=path + "TimeQualityViolinPlot.html",
                     title="Violin plot of quality over time")

    fig = go.Figure()

    fig.add_trace(go.Violin(y=dfs["quals"],
                            x=dfs["timebin"],
                            points=False, spanmode="hard",
                            line_color='black', line_width=1.5,
                            fillcolor=color, opacity=0.8))

    fig.update_layout(xaxis_title='Interval (hours)',
                      yaxis_title='Basecall quality',
                      title=title or time_qual.title,
                      title_x=0.5)

    fig.update_xaxes(tickangle=45)

    time_qual.fig = fig
    time_qual.html = time_qual.fig.to_html(full_html=False, include_plotlyjs='cdn')
    time_qual.save(figformat)

    return time_qual


def sequencing_speed_over_time(dfs, path, title, figformat, color="#4CB391"):
    time_duration = Plot(path=path + "TimeSequencingSpeed_ViolinPlot.html",
                         title="Violin plot of sequencing speed over time")

    mask = dfs['duration'] != 0

    fig = go.Figure()

    fig.add_trace(
        go.Violin(x=dfs.loc[mask, "timebin"],
                  y=dfs.loc[mask, "lengths"] / dfs.loc[mask, "duration"],
                  points=False, spanmode="hard",
                  line_color='black', line_width=1.5,
                  fillcolor=color, opacity=0.8))

    fig.update_layout(xaxis_title='Interval (hours)',
                      yaxis_title='Sequencing speed (nucleotides/second)',
                      title=title or time_duration.title,
                      title_x=0.5)

    fig.update_xaxes(tickangle=45)

    time_duration.fig = fig
    time_duration.html = time_duration.fig.to_html(full_html=False, include_plotlyjs='cdn')
    time_duration.save(figformat)

    return time_duration


def add_time_bins(dfs, bin_length=3):
    maxtime = dfs["start_time"].max().total_seconds()
    labels = [str(i) + "-" + str(i + bin_length)
              for i in range(0, 168, bin_length) if not i >= (maxtime / 3600)]
    return pd.cut(x=dfs["start_time"],
                  bins=ceil((maxtime / 3600) / bin_length),
                  labels=labels)


def plot_over_time(dfs, path, title, figformat, color="#4CB391"):
    num_reads = Plot(path=path + "NumberOfReads_Over_Time.html",
                     title="Number of reads over time")
    s = dfs.loc[:, "lengths"].resample('10T').count()

    fig = px.scatter(
        data_frame=None,
        x=s.index.total_seconds() / 3600,
        y=s)
    fig.update_traces(marker=dict(color=color))

    fig.update_layout(xaxis_title='Run time (hours)',
                      yaxis_title='Number of reads per 10 minutes',
                      title=title or num_reads.title,
                      title_x=0.5)

    num_reads.fig = fig
    num_reads.html = num_reads.fig.to_html(full_html=False, include_plotlyjs='cdn')
    num_reads.save(figformat)

    plots = [num_reads]

    if "channelIDs" in dfs:
        pores_over_time = Plot(path=path + "ActivePores_Over_Time.html",
                               title="Number of active pores over time")
        s = dfs.loc[:, "channelIDs"].resample('10T').nunique()

        fig = px.scatter(
            data_frame=None,
            x=s.index.total_seconds() / 3600,
            y=s)
        fig.update_traces(marker=dict(color=color))

        fig.update_layout(xaxis_title='Run time (hours)',
                          yaxis_title='Active pores per 10 minutes',
                          title=title or pores_over_time.title,
                          title_x=0.5)

        pores_over_time.fig = fig
        pores_over_time.html = pores_over_time.fig.to_html(full_html=False, include_plotlyjs='cdn')
        pores_over_time.save(figformat)

        plots.append(pores_over_time)
    return plots


def cumulative_yield(dfs, path, title, color, figformat):
    cum_yield_gb = Plot(path=path + "CumulativeYieldPlot_Gigabases.html",
                        title="Cumulative yield")

    s = dfs.loc[:, "lengths"].cumsum().resample('10T').max() / 1e9

    fig = px.scatter(
        x=s.index.total_seconds() / 3600,
        y=s)
    fig.update_traces(marker=dict(color=color))

    fig.update_layout(xaxis_title='Run time (hours)',
                      yaxis_title='Cumulative yield in gigabase',
                      title=title or cum_yield_gb.title,
                      title_x=0.5)

    cum_yield_gb.fig = fig
    cum_yield_gb.html = cum_yield_gb.fig.to_html(full_html=False, include_plotlyjs='cdn')
    cum_yield_gb.save(figformat)

    cum_yield_reads = Plot(path=path + "CumulativeYieldPlot_NumberOfReads.html",
                           title="Cumulative yield")

    s = dfs.loc[:, "lengths"].resample('10T').count().cumsum()

    fig = px.scatter(
        x=s.index.total_seconds() / 3600,
        y=s)
    fig.update_traces(marker=dict(color=color))

    fig.update_layout(xaxis_title='Run time (hours)',
                      yaxis_title='Cumulative yield in number of reads',
                      title=title or cum_yield_gb.title,
                      title_x=0.5)

    cum_yield_reads.fig = fig
    cum_yield_reads.html = cum_yield_reads.fig.to_html(full_html=False, include_plotlyjs='cdn')
    cum_yield_reads.save(figformat)

    return [cum_yield_gb, cum_yield_reads]
