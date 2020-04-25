import numpy as np
import logging
from nanoplotter.plot import Plot
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


class Layout(object):
    def __init__(self, structure, template, xticks, yticks, flowcell):
        self.structure = structure
        self.template = template
        self.xticks = xticks
        self.yticks = yticks
        self.flowcell = flowcell


def make_layout(maxval):
    """Make the physical layout of the MinION flowcell.
    based on https://bioinformatics.stackexchange.com/a/749/681
    returned as a numpy array
    """
    if maxval < 127:  # Flongle up to channel 126
        return Layout(
            structure=np.array([
                list(range(1, 13)) + [0],
                list(range(13, 25)) + [0],
                range(25, 38),
                range(38, 51),
                range(51, 64),
                range(64, 77),
                range(77, 90),
                range(90, 103),
                list(range(103, 115)) + [0],
                list(range(115, 127)) + [0]
            ]),
            template=np.zeros((10, 13)),
            xticks=range(1, 14),
            yticks=range(1, 11),
            flowcell='Flongle')
    elif maxval < 513:  # MinION up to channel 512
        layoutlist = []
        for i, j in zip(
                [33, 481, 417, 353, 289, 225, 161, 97],
                [8, 456, 392, 328, 264, 200, 136, 72]):
            for n in range(4):
                layoutlist.append(list(range(i + n * 8, (i + n * 8) + 8, 1))
                                  + list(range(j + n * 8, (j + n * 8) - 8, -1)))
        return Layout(
            structure=np.array(layoutlist).transpose(),
            template=np.zeros((16, 32)),
            xticks=range(1, 33),
            yticks=range(1, 17),
            flowcell='MinION')
    else:  # Assuming PromethION
        return Layout(
            structure=np.concatenate([np.array([list(range(10 * i + 1, i * 10 + 11))
                                                for i in range(25)]) + j
                                      for j in range(0, 3000, 250)],
                                     axis=1),
            template=np.zeros((25, 120)),
            xticks=range(1, 121),
            yticks=range(1, 26),
            flowcell='PromethION')


def spatial_heatmap(array, path, title=None, color="Greens", figformat="png"):
    """Taking channel information and creating post run channel activity plots."""
    logging.info("Nanoplotter: Creating heatmap of reads per channel using {} reads."
                 .format(array.size))
    activity_map = Plot(
        path=path + "." + figformat,
        title="Number of reads generated per channel")
    layout = make_layout(maxval=np.amax(array))
    valueCounts = pd.value_counts(pd.Series(array))
    for entry in valueCounts.keys():
        layout.template[np.where(layout.structure == entry)] = valueCounts[entry]
    plt.figure()
    ax = sns.heatmap(
        data=pd.DataFrame(layout.template, index=layout.yticks, columns=layout.xticks),
        xticklabels="auto",
        yticklabels="auto",
        square=True,
        cbar_kws={"orientation": "horizontal"},
        cmap=color,
        linewidths=0.20)
    ax.set_title(title or activity_map.title)
    activity_map.fig = ax.get_figure()
    activity_map.save(format=figformat)
    plt.close("all")
    return [activity_map]
