import pandas as pd
import numpy as np


class BarcodeTitle(object):
    """Bit of a dummy class to add barcode titles to the report"""

    def __init__(self, title):
        self.title = title.upper()

    def encode(self):
        return ""


def chunks(values, chunks):
    if values:
        chunksize = int(len(values) / chunks)
        return([' '.join(values[i:i + chunksize]) for i in range(0, len(values), chunksize)])
    else:
        return [" "] * chunks


def html_stats(settings):
    statsfile = settings["statsfile"]
    filtered = settings["filtered"]
    as_tsv = settings['tsv_stats']

    stats_html = []
    stats_html.append('<div class="panel panelM"> <h1>NanoPlot report</h1>')
    if filtered:
        stats_html.append('<h2 id="stats0">Summary statistics prior to filtering</h2>')
        if as_tsv:
            stats_html.append(statsfile[0].to_html())
            stats_html.append('<h2 id="stats1">Summary statistics after filtering</h2>')
            stats_html.append(statsfile[1].to_html())
        else:
            stats_html.append(stats2html(statsfile[0]))
            stats_html.append('<h2 id="stats1">Summary statistics after filtering</h2>')
            stats_html.append(stats2html(statsfile[1]))
    else:
        stats_html.append('<h2 id="stats0">Summary statistics</h2>')
        if as_tsv:
            stats_html.append(statsfile[0].to_html())
        else:
            stats_html.append(stats2html(statsfile[0]))
    return '\n'.join(stats_html)


def stats2html(statsf):
    df = pd.read_csv(statsf, sep=':', header=None, names=['feature', 'value'])
    values = df["value"].str.strip().str.replace('\t', ' ').str.split().replace(np.nan, '')
    num = len(values[0]) or 1
    v = [chunks(i, num) for i in values]
    return pd.DataFrame(v, index=df["feature"]).to_html(header=False)


def html_toc(plots, filtered=False):
    toc = []
    toc.append('<div class="panel panelC">')
    if filtered:
        toc.append(
            '<p><strong><a href="#stats0">Summary Statistics prior to filtering</a></strong></p>')
        toc.append(
            '<p><strong><a href="#stats1">Summary Statistics after filtering</a></strong></p>')
    else:
        toc.append(
            '<p><strong><a href="#stats0">Summary Statistics</a></strong></p>')
    toc.append('<p><strong><a href="#plots">Plots</a></strong></p>')
    toc.extend(['<p style="margin-left:20px"><a href="#'
                + p.title.replace(' ', '_') + '">' + p.title + '</a></p>' for p in plots])
    toc.append('</div>')
    return '\n'.join(toc)


def html_plots(plots):
    html_plots = []
    html_plots.append('<h2 id="plots">Plots</h2>')
    for plot in plots:
        html_plots.append('\n<h3 id="' + plot.title.replace(' ', '_') + '">'
                          + plot.title + '</h3>\n' + plot.encode())
        html_plots.append('\n<br>\n<br>\n<br>\n<br>')
    return '\n'.join(html_plots)


def run_info(settings):
    html_info = []
    html_info.append('<h4>Run Info</h4>\n')
    html_info.append('<h4>Data source:</h5>\n')
    for k in ["fastq", "fasta", "fastq_rich", "fastq_minimal", "summary",
              "bam", "ubam", "cram", "pickle", "feather"]:
        html_info.append(f"{k}:\t{settings[k]}<br>")
    html_info.append('<h4>Filtering parameters:</h5>\n')
    for k in ['maxlength', 'minlength', 'drop_outliers', 'downsample', 'loglength',
              'percentqual', 'alength', 'minqual', 'runtime_until', 'no_supplementary']:
        html_info.append(f"{k}:\t{settings[k]}<br>")
    # html_info.append('</p>')
    return '\n'.join(html_info)


html_head = """<!DOCTYPE html>
<html>
    <head>
    <meta charset="UTF-8">
        <style>
        table, th, td {
            text-align: left;
            padding: 2px;
            /* border: 1px solid black;
            border-collapse: collapse; */
        }
        h2 {
            line-height: 0pt;
        }
        .panel {
            display: inline-block;
            background: #ffffff;
            min-height: 100px;
            box-shadow:0px 0px 5px 5px #C9C9C9;
            -webkit-box-shadow:2px 2px 5px 5x #C9C9C9;
            -moz-box-shadow:2px 2px 5px 5px #C9C9C9;
            margin: 10px;
            padding: 10px;
        }
        .panelC {
            float: left
        }
        .panelM {
            float: left
        }
        </style>
        <title>NanoPlot Report</title>
    </head>"""
