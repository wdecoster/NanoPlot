# wdecoster
"""
This module provides a few simple math and statistics functions
for other scripts processing Oxford Nanopore sequencing data




# FUNCTIONS
* Calculate read N50 from a set of lengths
get_N50(readlenghts)
* Remove extreme length outliers from a dataset
remove_length_outliers(dataframe, columname)
* Calculate the average Phred quality of a read
ave_qual(qualscores)
* Write out the statistics report after calling readstats function
write_stats(dataframe, outputname)
* Compute a number of statistics, return a dictionary
calc_read_stats(dataframe)
"""

import numpy as np
import sys
from math import log


class Stats(object):
    def __init__(self, df):
        if len(df) < 5:
            sys.stderr.write("\n\nWARNING: less than 5 reads in the dataset!\n")
            sys.stderr.write("WARNING: some stats might be unexpected or missing\n")
            sys.stderr.write("WARNING: or a crash might happen, who knows\n")
            sys.stderr.write(
                "WARNING: this code is not intended for such small datasets\n\n\n"
            )
        self.number_of_reads = len(df)
        self.number_of_bases = np.sum(df["lengths"])
        self._with_readIDs = "readIDs" in df
        if "aligned_lengths" in df:
            self.number_of_bases_aligned = np.sum(df["aligned_lengths"])
            self.fraction_bases_aligned = (
                self.number_of_bases_aligned / self.number_of_bases
            )
        self.median_read_length = np.median(df["lengths"])
        self.mean_read_length = np.mean(df["lengths"])
        self.read_length_stdev = np.std(df["lengths"])
        self.n50 = get_N50(np.sort(df["lengths"]))
        if "percentIdentity" in df:
            self.average_identity = np.mean(df["percentIdentity"])
            self.median_identity = np.median(df["percentIdentity"])
        if "channelIDs" in df:
            self.active_channels = np.unique(df["channelIDs"]).size
        if "quals" in df:
            self._qualgroups = [
                10,
                15,
                20,
                25,
                30,
            ]  # needs 5 elements in current implementation
            self.mean_qual = ave_qual_floats(df["quals"])
            self.median_qual = np.median(df["quals"])
            self._top5_lengths = get_top_5(
                df=df, col="lengths", values=["lengths", "quals"]
            )
            self._top5_quals = get_top_5(
                df=df, col="quals", values=["quals", "lengths"]
            )
            self._reads_above_qual = [reads_above_qual(df, q) for q in self._qualgroups]
        else:
            self._top5_lengths = get_top_5(
                df=df, col="lengths", values=["lengths"], fill="quals"
            )

    def long_features_as_string(self):
        """formatting long features to a string to print for legacy stats output"""
        self.top5_lengths = self.long_feature_as_string_top5(self._top5_lengths)
        self.top5_quals = self.long_feature_as_string_top5(self._top5_quals)
        self.reads_above_qual = self.long_feature_as_string_above_qual(
            self._reads_above_qual
        )

    def long_feature_as_string_top5(self, field):
        """for legacy stats output"""
        if self._with_readIDs:
            return [
                str(round(i, ndigits=1))
                + " ("
                + str(round(j, ndigits=1))
                + "; "
                + k
                + ")"
                for i, j, k in field
            ]
        else:
            return [
                str(round(i, ndigits=1)) + " (" + str(round(j, ndigits=1)) + ")"
                for i, j in field
            ]

    def long_feature_as_string_above_qual(self, field):
        """for legacy stats output"""
        return [self.format_above_qual_line(entry) for entry in field]

    def format_above_qual_line(self, entry):
        """for legacy stats output"""
        numberAboveQ, megAboveQ = entry
        return "{} ({}%) {}Mb".format(
            numberAboveQ,
            round(100 * (numberAboveQ / self.number_of_reads), ndigits=1),
            round(megAboveQ, ndigits=1),
        )

    def to_dict(self):
        """for tsv stats output"""
        statdict = self.__dict__
        for key, value in statdict.items():
            if not key.startswith("_"):
                if not isinstance(value, int):
                    statdict[key] = "{:.1f}".format(value)
        self.unwind_long_features_top5(
            feature="_top5_lengths", name="longest_read_(with_Q)"
        )
        self.unwind_long_features_top5(
            feature="_top5_quals", name="highest_Q_read_(with_length)"
        )
        self.unwind_long_features_above_qual(feature="_reads_above_qual", name="Reads")
        return {k: v for k, v in statdict.items() if not k.startswith("_")}

    def unwind_long_features_top5(self, feature, name):
        """for tsv stats output"""
        if feature not in self.__dict__:
            return
        for entry, label in zip(self.__dict__[feature], range(1, 6)):
            self.__dict__[name + ":" + str(label)] = "{} ({})".format(
                round(entry[0], ndigits=1), round(entry[1], ndigits=1)
            )

    def unwind_long_features_above_qual(self, feature, name):
        """for tsv stats output"""
        if feature not in self.__dict__:
            return
        for entry, label in zip(
            self.__dict__[feature], [">Q{}:".format(q) for q in self._qualgroups]
        ):
            numberAboveQ, megAboveQ = entry
            percentage = 100 * (numberAboveQ / float(self.number_of_reads))
            self.__dict__[name + " " + label] = "{} ({}%) {}Mb".format(
                numberAboveQ, round(percentage, ndigits=1), round(megAboveQ, ndigits=1)
            )


def get_N50(readlengths):
    """Calculate read length N50.

    Based on https://github.com/PapenfussLab/Mungo/blob/master/bin/fasta_stats.py
    """
    return readlengths[
        np.where(np.cumsum(readlengths) >= 0.5 * np.sum(readlengths))[0][0]
    ]


def errs_tab(n):
    """Generate list of error rates for qualities less than equal than n."""
    return [10 ** (q / -10) for q in range(n + 1)]


def ave_qual(quals, qround=False, tab=errs_tab(128)):
    """Calculate average basecall quality of a read.

    Receive the integer quality scores of a read and return the average quality for that read
    First convert Phred scores to probabilities,
    calculate average error probability
    convert average back to Phred scale
    """
    if quals:
        mq = -10 * log(sum([tab[q] for q in quals]) / len(quals), 10)
        if qround:
            return round(mq)
        else:
            return mq
    else:
        return None

def ave_qual_floats(quals):
    """
    This function is to create the average quality across reads, where the input is a float
    """
    convert_to_probs = lambda q: 10 ** (-q / 10)
    vfunc = np.vectorize(convert_to_probs)
    probs = vfunc(quals)
    return -10 * log(probs.sum() / len(probs), 10)

def get_top_5(df, col, values, fill=False):
    if "readIDs" in df:
        values.append("readIDs")
    if fill:
        return (
            df.sort_values(col, ascending=False)
            .head(5)[values]
            .assign(fill=[0] * 5)
            .reset_index(drop=True)
            .itertuples(index=False, name=None)
        )
    else:
        return (
            df.sort_values(col, ascending=False)
            .head(5)[values]
            .reset_index(drop=True)
            .itertuples(index=False, name=None)
        )


def reads_above_qual(df, qual):
    numberAboveQ = np.sum(df["quals"] > qual)
    megAboveQ = np.sum(df.loc[df["quals"] > qual, "lengths"]) / 1e6
    return numberAboveQ, megAboveQ


def write_stats(datadfs, outputfile, names=[], as_tsv=False):
    """Call calculation functions and write stats file.

    This function takes a list of DataFrames,
    and will create a column for each in the tab separated output.
    """
    if outputfile == "stdout":
        output = sys.stdout
    else:
        output = open(outputfile, "wt")

    stats = [Stats(df) for df in datadfs]

    if as_tsv:
        import pandas as pd

        df = pd.DataFrame([s.to_dict() for s in stats]).transpose()
        df.index.name = "Metrics"
        if names:
            df.columns = names
        else:
            df.columns = ["dataset"]
        output.write(df.to_csv(sep="\t"))
        return df
    else:
        write_stats_legacy(stats, names, output, datadfs)


def write_stats_legacy(stats, names, output, datadfs):
    """
    Legacy method to write out stats.
    Will add padding to pretty print the table, and contain section headers
    """
    features = {
        "Number of reads": "number_of_reads",
        "Total bases": "number_of_bases",
        "Total bases aligned": "number_of_bases_aligned",
        "Fraction of bases aligned": "fraction_bases_aligned",
        "Median read length": "median_read_length",
        "Mean read length": "mean_read_length",
        "STDEV read length": "read_length_stdev",
        "Read length N50": "n50",
        "Average percent identity": "average_identity",
        "Median percent identity": "median_identity",
        "Active channels": "active_channels",
        "Mean read quality": "mean_qual",
        "Median read quality": "median_qual",
    }
    max_len = max([len(k) for k in features.keys()])
    try:
        max_num = (
            max(
                max([len(str(s.number_of_bases)) for s in stats]),
                max([len(str(n)) for n in names]),
            )
            + 6
        )
    except ValueError:
        max_num = max([len(str(s.number_of_bases)) for s in stats]) + 6
    output.write(
        "{:<{}}{}\n".format(
            "General summary:",
            max_len,
            " ".join(["{:>{}}".format(n, max_num) for n in names]),
        )
    )
    for f in sorted(features.keys()):
        try:
            output.write(
                "{f:{pad}}{v}\n".format(
                    f=f + ":",
                    pad=max_len,
                    v=feature_list(stats, features[f], padding=max_num),
                )
            )
        except KeyError:
            pass
    if all(["quals" in df for df in datadfs]):
        for s in stats:
            s.long_features_as_string()
        long_features = {
            "Top 5 longest reads and their mean basecall quality score": [
                "top5_lengths",
                range(1, 6),
            ],
            "Top 5 highest mean basecall quality scores and their read lengths": [
                "top5_quals",
                range(1, 6),
            ],
            "Number, percentage and megabases of reads above quality cutoffs": [
                "reads_above_qual",
                [">Q" + str(q) for q in stats[0]._qualgroups],
            ],
        }
        for lf in sorted(long_features.keys()):
            output.write(lf + "\n")
            for index in range(5):
                output.write(
                    "{}:\t{}\n".format(
                        long_features[lf][1][index],
                        feature_list(
                            stats=stats, feature=long_features[lf][0], index=index
                        ),
                    )
                )


def feature_list(stats, feature, index=None, padding=15):
    if index is None:
        return " ".join(
            ["{:>{},.1f}".format(s.__dict__[feature], padding) for s in stats]
        )
    else:
        return "\t".join(
            [
                str(s.__dict__[feature][index])
                if len(s.__dict__[feature]) > index
                else "NA"
                for s in stats
            ]
        )
