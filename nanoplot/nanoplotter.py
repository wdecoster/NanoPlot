# wdecoster

from __future__ import division
import time
import logging
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

def scatter(x, y, names, path, stat=None, log=False, minvalx=0, minvaly=0):
	'''
	Plotting function
	Create three types of joint plots of length vs quality, containing marginal summaries
	-A scatter plot with histograms on axes
	-A hexagonal binned plot with histograms on axes
	-A kernel density plot with density curves on axes, subsampled to 10000 reads if required
	'''
	logging.info("Creating {} vs {} plots using statistics from {} reads.".format(names[0], names[1], x.size))
	sns.set(style="ticks")
	maxvalx = np.amax(x)
	maxvaly = np.amax(y)
	plot = sns.jointplot(
		x=x,
		y=y,
		kind="hex",
		color="#4CB391",
		stat_func=stat,
		space=0,
		xlim=(minvalx, maxvalx),
		ylim=(minvaly, maxvaly),
		size=10)
	plot.set_axis_labels(names[0], names[1])
	if log:
		plot.ax_joint.set_xticklabels(10**plot.ax_joint.get_xticks().astype(int))
	plot.savefig(path + "_hex.png", format='png', dpi=1000)
	sns.set(style="darkgrid")
	plot = sns.jointplot(
		x=x,
		y=y,
		kind="scatter",
		color="#4CB391",
		stat_func=stat,
		xlim=(minvalx, maxvalx),
		ylim=(minvaly, maxvaly),
		space=0,
		size=10,
		joint_kws={"s": 1})
	plot.set_axis_labels(names[0], names[1])
	if log:
		plot.ax_joint.set_xticklabels(10**plot.ax_joint.get_xticks().astype(int))
	plot.savefig(path + "_scatter.png", format='png', dpi=1000)

	plot = sns.jointplot(
		x=x,
		y=y,
		kind="kde",
		clip=((0, np.Inf), (0, np.Inf)),
		xlim=(minvalx, maxvalx),
		ylim=(minvaly, maxvaly),
		space=0,
		color="#4CB391",
		stat_func=stat,
		shade_lowest=False,
		size=10)
	plot.set_axis_labels(names[0], names[1])
	if log:
		plot.ax_joint.set_xticklabels(10**plot.ax_joint.get_xticks().astype(int))
	plot.savefig(path + "_kde.png", format='png', dpi=1000)
	plt.close("all")


def timePlots(df, path):
	'''
	Plotting function
	Making plots of time vs read length, time vs quality and cumulative yield
	'''
	dfs_sparse = df.sample(min(2000, len(df.index))).sort_values("start_time")
	dfs_sparse["time_h"] = dfs_sparse["start_time"].astype('timedelta64[h]')
	g = sns.JointGrid(
		x='time_h',
		y='quals',
		data=dfs_sparse,
		space=0,
		size=10,
		xlim=(0, dfs_sparse.time_h.max()))
	g.plot_joint(plt.scatter, color="#4CB391")
	g.ax_marg_y.hist(dfs_sparse['quals'].dropna(), orientation="horizontal", color="#4CB391")
	g.set_axis_labels('Run tim (hours)', 'Median average basecall quality')
	g.savefig(path + "TimeQualityScatterPlot.png", format='png', dpi=1000)

	g = sns.JointGrid(
		x='time_h',
		y="lengths",
		data=dfs_sparse,
		space=0,
		size=10,
		xlim=(0, dfs_sparse.time_h.max()))
	g.plot_joint(plt.scatter, color="#4CB391")
	g.ax_marg_y.hist(dfs_sparse["lengths"].dropna(), orientation="horizontal", color="#4CB391")
	g.set_axis_labels('Run tim (hours)', 'Median read length')
	g.savefig(path + "TimeLengthScatterPlot.png", format='png', dpi=1000)

	dfs_sparse["cumyield_gb"] = dfs_sparse["lengths"].cumsum() / 10**9
	g = sns.JointGrid(
		x='time_h',
		y="cumyield_gb",
		data=dfs_sparse,
		space=0,
		size=10,
		xlim=(0, dfs_sparse.time_h.max()),
		ylim=(0, dfs_sparse.tail(1)["cumyield_gb"].item()))
	g.plot_joint(plt.scatter, color="#4CB391")
	g.set_axis_labels('Run tim (hours)', 'Cumulative yield in gigabase')
	g.savefig(path + "CumulativeYieldPlot.png", format='png', dpi=1000)
	plt.close("all")


def lengthPlots(array, name, path, n50, log=False):
	'''
	Plotting function
	Create density plot and histogram based on a numpy array (read lengths or transformed read lengths)
	'''
	import math
	maxvalx = np.amax(array)
	if log:
		logging.info("Creating length plots for {} from {} reads with read length N50 of {}.".format(name, array.size, n50))
	else:
		logging.info("Creating length plots for {} from {} reads with read length N50 of {}.".format(name, array.size, n50))

	ax = sns.distplot(
		a=array,
		kde=True,
		hist=False,
		color="#4CB391",
		kde_kws={"label": name, "clip": (0, maxvalx)})
	if log:
		ticks = [10**i for i in range(10) if not 10**i > 10**math.ceil(math.log(10**maxvalx,10))]
		ax.set(
			xticks=np.log10(ticks),
			xticklabels=ticks)
	fig = ax.get_figure()
	fig.savefig(path + "DensityCurve" + name.replace(' ', '') + ".png", format='png', dpi=1000)
	plt.close("all")

	ax = sns.distplot(
		a=array,
		kde=False,
		hist=True,
		color="#4CB391")

	if log:
		ax.set(
			xticks=np.log10(ticks),
			xticklabels=ticks)
		plt.axvline(np.log10(n50))
		plt.annotate('N50', xy=(np.log10(n50), np.amax([h.get_height() for h in ax.patches])), size=8)
	else:
		plt.axvline(n50)
		plt.annotate('N50', xy=(n50, np.amax([h.get_height() for h in ax.patches])), size=8)
	ax.set(xlabel='Read length', ylabel='Number of reads')
	fig = ax.get_figure()
	fig.savefig(path + "Histogram" + name.replace(' ', '') + ".png", format='png', dpi=1000)
	plt.close("all")


def spatialHeatmap(array, title, path, colour):
	'''
	Plotting function
	Taking channel information and creating post run channel activity plots
	'''
	logging.info("Creating activity maps for {} using statistics from {} reads.".format(title.lower(), array.size))
	layout = np.array([
		range(125, 0, -4), range(126, 1, -4), range(127, 2, -4), range(128, 3, -4),
		range(253, 128, -4), range(254, 129, -4), range(255, 130, -4), range(256, 131, -4),
		range(381, 256, -4), range(382, 257, -4), range(383, 258, -4), range(384, 259, -4),
		range(509, 384, -4), range(510, 385, -4), range(511, 386, -4), range(512, 387, -4)])
	activityData = np.zeros((16, 32))
	valueCounts = pd.value_counts(pd.Series(array))
	for entry in valueCounts.keys():
		activityData[np.where(layout == entry)] = valueCounts[entry]
	plt.figure()
	ax = sns.heatmap(
		data=activityData,
		xticklabels=range(1, 33),
		yticklabels=range(1, 17),
		square=True,
		cbar_kws={"orientation": "horizontal"},
		cmap=colour,
		linewidths=0.20)
	ax.set_title(title)
	fig = ax.get_figure()
	fig.savefig(path + ".png", format='png', dpi=1000)
	plt.close("all")
