import logging
import numpy as np
from datetime import timedelta
import sys


def flag_length_outliers(df, columnname):
    """Return index of records with length-outliers above 3 standard deviations from the median."""
    return df[columnname] > (np.median(df[columnname]) + 3 * np.std(df[columnname]))


def phred_to_percent(phred):
    return 100 * (1 - 10 ** (phred / -10))


def non_filtered_reads(df):
    return len(df[df["length_filter"]])


def filter_and_transform_data(df, settings):
    '''
    Perform filtering on the data based on arguments set on commandline
    - use aligned length or sequenced length (bam mode only)
    - hide outliers from length plots*
    - hide reads longer than maxlength or shorter than minlength from length plots*
    - filter reads with a quality below minqual
    - use log10 scaled reads rather than normal
    - use empirical percent accuracy rather than phred score quality
    - downsample reads to args.downsample

    - always: drop reads which are basecaller artefacts
              judged by length below 20 and quality above 30

    * using a boolean column length_filter
    '''
    df["length_filter"] = True
    settings["filtered"] = False

    if settings.get("alength") and settings.get("bam"):
        settings["lengths_pointer"] = "aligned_lengths"
        logging.info("Using aligned read lengths for plotting.")
    else:
        settings["lengths_pointer"] = "lengths"
        logging.info("Using sequenced read lengths for plotting.")

    if settings.get("drop_outliers"):
        num_reads_prior = non_filtered_reads(df)
        df.loc[flag_length_outliers(df, settings["lengths_pointer"]), "length_filter"] = False
        num_reads_post = non_filtered_reads(df)
        logging.info("Hidding {} length outliers in length plots.".format(
            str(num_reads_prior - num_reads_post)))

    if settings.get("maxlength"):
        num_reads_prior = non_filtered_reads(df)
        df.loc[df[settings["lengths_pointer"]] > settings["maxlength"], "length_filter"] = False
        num_reads_post = non_filtered_reads(df)
        logging.info("Hidding {} reads longer than {}bp in length plots.".format(
            str(num_reads_prior - num_reads_post),
            str(settings["maxlength"])))

    if settings.get("minlength"):
        num_reads_prior = non_filtered_reads(df)
        df.loc[df[settings["lengths_pointer"]] < settings["minlength"], "length_filter"] = False
        num_reads_post = non_filtered_reads(df)
        logging.info("Hidding {} reads shorter than {}bp in length plots.".format(
            str(num_reads_prior - num_reads_post),
            str(settings["minlength"])))

    if settings.get("minqual"):
        if "quals" in df:
            num_reads_prior = non_filtered_reads(df)
            df = df.loc[df["quals"] > settings["minqual"]].copy()
            num_reads_post = non_filtered_reads(df)
            logging.info("Removing {} reads with quality below Q{}.".format(
                str(num_reads_prior - num_reads_post),
                str(settings["minqual"])))
            settings["filtered"] = True
        else:
            sys.stderr.write("--minqual is ignored since no quality information in the data.")
            logging.info("--minqual is ignored since no quality information in the data.")

    if settings.get("loglength"):
        df["log_" + settings["lengths_pointer"]] = np.log10(df[settings["lengths_pointer"]])
        settings["lengths_pointer"] = "log_" + settings["lengths_pointer"]
        logging.info("Using log10 scaled read lengths.")
        settings["logBool"] = True
    else:
        settings["logBool"] = False

    if settings.get("runtime_until"):
        if "start_time" in df:
            num_reads_prior = non_filtered_reads(df)
            df = df[df.start_time < timedelta(hours=settings["runtime_until"])]
            num_reads_post = non_filtered_reads(df)
            logging.info("Removing {} reads generated after {} hours in the run.".format(
                str(num_reads_prior - num_reads_post),
                str(settings["runtime_until"])))
            settings["filtered"] = True
        else:
            sys.stderr.write("--runtime_until is ignored since no time information in the data.")
            logging.info("--runtime_until is ignored since no time information in the data.")

    if "quals" in df:
        num_reads_prior = len(df)
        df = df.loc[-((df["lengths"] < 20) & (df["quals"] > 30))].copy()
        num_reads_post = len(df)
        if num_reads_prior - num_reads_post > 0:
            logging.info(
                "Removed {} artefactual reads with very short length and very high quality."
                .format(num_reads_prior - num_reads_post))
            settings["filtered"] = True

    if settings.get("downsample"):
        new_size = min(settings["downsample"], len(df))
        logging.info("Downsampling the dataset from {} to {} reads".format(
            len(df), new_size))
        df = df.sample(new_size)
        settings["filtered"] = True

    if settings.get("percentqual"):
        df["quals"] = df["quals"].apply(phred_to_percent)
        logging.info("Converting quality scores to theoretical percent identities.")

    return(df, settings)
