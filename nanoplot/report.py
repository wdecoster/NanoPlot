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
        return ([' '.join(values[i:i + chunksize]) for i in range(0, len(values), chunksize)])
    else:
        return [" "] * chunks


def html_stats(settings):
    statsfile = settings["statsfile"]
    filtered = settings["filtered"]
    as_tsv = settings['tsv_stats']

    stats_html = []
    stats_html.append('<main class="grid-main"><h2>NanoPlot reports</h2>')
    if filtered:
        stats_html.append('<h3 id="stats0">Summary statistics prior to filtering</h3>')
        if as_tsv:
            stats_html.append(statsfile[0].to_html())
            stats_html.append('<h3 id="stats1">Summary statistics after filtering</h3>')
            stats_html.append(statsfile[1].to_html())
        else:
            stats_html.append(stats2html(statsfile[0]))
            stats_html.append('<h3 id="stats1">Summary statistics after filtering</h3>')
            stats_html.append(stats2html(statsfile[1]))
    else:
        stats_html.append('<h3 id="stats0">Summary statistics</h3>')
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
    df = pd.DataFrame(v, index=df["feature"])
    df.columns.name = None
    df.index.name = None
    return df.to_html(header=False)


def html_toc(plots, filtered=False):
    toc = []
    toc.append('<h1 class="hiddentitle">NanoPlot statistics report</h1>')
    toc.append('<header class="grid-header"><nav><h2 class="hiddentitle">Menu</h2><ul>')
    if filtered:
        toc.append(
            '<li><a href="#stats0">Summary Statistics prior to filtering</a></li>')
        toc.append(
            '<li><a href="#stats1">Summary Statistics after filtering</a></li>')
    else:
        toc.append('<li><a href="#stats0">Summary Statistics</a></li>')

    toc.append('<li class="submenu"><a href="#plots" class="submenubtn">Plots</a>')
    toc.append('<ul class="submenu-items">')
    toc.extend(['<li><a href="#'
                + p.title.replace(' ', '_') + '">' + p.title + '</a></li>' for p in plots])
    toc.append('</ul>')
    toc.append('</li>')
    toc.append(
        '<li class="issue-btn"><a href="https://github.com/wdecoster/NanoPlot/issues" target="_blank"  class="reporting">Report issue on Github</a></li>')
    toc.append('</ul></nav></header>')
    return '\n'.join(toc)


def html_plots(plots):
    html_plots = []
    html_plots.append('<h3 id="plots">Plots</h3>')
    for plot in plots:
        html_plots.append('<button class="collapsible">' + plot.title + '</button>')
        html_plots.append('<section class="collapsible-content"><h4 class="hiddentitle" id="' +
                          plot.title.replace(' ', '_') + '">' + plot.title + '</h4>')
        html_plots.append(plot.encode())
        html_plots.append('</section>')

    html_plots.append(
        '<script>var coll = document.getElementsByClassName("collapsible");var i;for (i = 0; i < coll.length; i++) {coll[i].addEventListener("click", function() {this.classList.toggle("active");var content = this.nextElementSibling;if (content.style.display === "none") {content.style.display = "block";} else {content.style.display = "none";}});}</script>')

    return '\n'.join(html_plots)


def run_info(settings):
    html_info = []
    html_info.append('<h5>Run Info</h5>\n')
    html_info.append('<h6>Data source:</h6>\n')
    for k in ["fastq", "fasta", "fastq_rich", "fastq_minimal", "summary",
              "bam", "ubam", "cram", "pickle", "feather"]:
        html_info.append(f"{k}:\t{settings[k]}<br>")
    html_info.append('<h6>Filtering parameters:</h6>\n')
    for k in ['maxlength', 'minlength', 'drop_outliers', 'downsample', 'loglength',
              'percentqual', 'alength', 'minqual', 'runtime_until', 'no_supplementary']:
        html_info.append(f"{k}:\t{settings[k]}<br>")
    # html_info.append('</p>')
    return '\n'.join(html_info)


html_head = """
<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<style>

body {margin:0}

.grid { /* grid definition for index page */
    display: grid;
    grid-template-areas:    'gheader'
                            'gmain';
    margin: 0;
}

.grid > .grid-header { /* definition of the header on index page and its position in the grid */
    grid-area: gheader;
}

.grid > .grid-main { /* definition of the main content on index page and its position in the grid */
    grid-area: gmain;
}

nav {
    text-align: center;
}

ul {
    border-bottom: 1px solid white;
    font-family: "Trebuchet MS", sans-serif;
    list-style-type: none; /* remove dot symbols from list */
    margin: 0;
    padding: 0;
    overflow: hidden; /* contains the overflow of the element if it goes 'out of bounds' */
    background-color: #001f3f;
    font-size: 1.6em;
}

ul > li > ul {
    font-size: 1em;
}

li {
    float: left; /* floats the list items to the left side of the page */
}

li a, .submenubutton {
    display: inline-block; /* display the list items inline block so the items are vertically displayed */
    color: white;
    text-align: center;
    padding: 14px 16px;
    text-decoration: none; /* removes the underline that comes with the a tag */
}

li a:hover, .submenu:hover .submenubutton { /* when you hover over a submenu item the bkgrnd color is gray */
    background-color: #39CCCC;
}

.submenu {
    display: inline-block; /* idem to above, list items are displayed underneath each other */
}

.submenu-items { /* hides the ul */
    display: none;
    position: absolute;
    background-color: #f9f9f9;
    min-width: 160px;
    z-index: 1;
}

.submenu-items li {
    display: block;
    float: none;
    overflow: hidden;
}

.submenu-items li a { /* styling of the links in the submenu */
    color: black;
    padding: 12px 16px;
    text-decoration: none;
    display: block;
    text-align: left;
}

.submenu-items a:hover {
    background-color: #f1f1f1;
}

.submenu:hover .submenu-items {
    display: block;
    float: bottom;
    overflow: hidden;
}

li {
  border-right: 1px solid #bbb;
}

.issue-btn {
  border-right: none;
  float: right;
}

.hiddentitle { /* hides titles that are not necessary for content, but are for outline */
  position: absolute;
  width: 1px;
  height: 1px;
  overflow: hidden;
  left: -10000px;
}

h2 { color: #111; font-family: 'Helvetica Neue', sans-serif; font-size: 60px; font-weight: bold; letter-spacing: -1px; line-height: 1; text-align: center; }

h3 { color: #111; font-family: 'Open Sans', sans-serif; font-size: 25px; font-weight: 300; line-height: 32px; text-align: center; padding-bottom: 0;}

h4 { color: #111; font-family: 'Helvetica Neue', sans-serif; font-size: 16px; font-weight: 150; margin: 0 0 0 0; text-align: left; padding:20px 0px 20px 0px;}

table {
  font-family: Arial, Helvetica, sans-serif;
  border-collapse: collapse;
  table-layout: auto;
  border-collapse: collapse;
  width: 100%;
}

table td, table th {
  border: 1px solid #ddd;
  padding: 8px;
}

table tr:nth-child(even){background-color: #f2f2f2;}

table tr:hover {background-color: #ddd;}

/* Style the button that is used to open and close the collapsible content */
.collapsible {
  background-color: #39CCCC;
  color: white;
  cursor: pointer;
  padding: 18px;
  width: 100%;
  border: none;
  text-align: left;
  outline: none;
  font-size: 15px;
}

/* Add a background color to the button if it is clicked on (add the .active class with JS), and when you move the mouse over it (hover) */
.active, .collapsible:hover {
    color:white;
  background-color: #001f3f;
}

/* Style the collapsible content. Note: hidden by default */
.collapsible-content {
  padding: 0 18px;
  display: block;
  overflow: hidden;
  background-color: #FFFFFF;
  text-align: center;
}

.collapsible:after {
  content: '-';
  font-size: 20px;
    font-weight: bold;
  float: right;
    color:white;
  margin-left: 5px;
}

.active:after {
  content: '+'; /* Unicode character for "minus" sign (-) */
      color: white;
}
</style>
<title>NanoPlot Report</title>
</head>
"""
