import logging
import nanomath
import numpy as np
from datetime import timedelta


def filter_and_transform_data(datadf, settings):
    '''
    Perform filtering on the data based on arguments set on commandline
    - use aligned length or sequenced length (bam mode only)
    - drop outliers
    - drop reads longer than args.maxlength
    - use log10 scaled reads
    - downsample reads to args.downsample
    Return an accurate prefix which is added to plotnames using this filtered data
    '''
    length_prefix_list = list()
    settings["filtered"] = False
    if settings.get("alength") and settings.get("bam"):
        settings["lengths_pointer"] = "aligned_lengths"
        length_prefix_list.append("Aligned_")
        logging.info("Using aligned read lengths for plotting.")
    else:
        settings["lengths_pointer"] = "lengths"
        logging.info("Using sequenced read lengths for plotting.")
    if settings.get("drop_outliers"):
        num_reads_prior = len(datadf)
        datadf = nanomath.remove_length_outliers(datadf, settings["lengths_pointer"])
        length_prefix_list.append("OutliersRemoved_")
        num_reads_post = len(datadf)
        logging.info("Removing {} length outliers for plotting.".format(
            str(num_reads_prior - num_reads_post)))
        settings["filtered"] = True
    if settings.get("maxlength"):
        num_reads_prior = len(datadf)
        datadf = datadf.loc[datadf[settings["lengths_pointer"]] < settings["maxlength"]].copy()
        length_prefix_list.append("MaxLength-" + str(settings["maxlength"]) + '_')
        num_reads_post = len(datadf)
        logging.info("Removed {} reads longer than {}bp.".format(
            str(num_reads_prior - num_reads_post),
            str(settings["maxlength"])))
        settings["filtered"] = True
    if settings.get("minlength"):
        num_reads_prior = len(datadf)
        datadf = datadf.loc[datadf[settings["lengths_pointer"]] > settings["minlength"]].copy()
        length_prefix_list.append("MinLength-" + str(settings["minlength"]) + '_')
        num_reads_post = len(datadf)
        logging.info("Removed {} reads shorter than {}bp.".format(
            str(num_reads_prior - num_reads_post),
            str(settings["minlength"])))
        settings["filtered"] = True
    if settings.get("minqual"):
        num_reads_prior = len(datadf)
        datadf = datadf.loc[datadf["quals"] > settings["minqual"]].copy()
        num_reads_post = len(datadf)
        logging.info("Removing {} reads with quality below Q{}.".format(
            str(num_reads_prior - num_reads_post),
            str(settings["minqual"])))
        settings["filtered"] = True
    if settings.get("loglength"):
        datadf["log_" + settings["lengths_pointer"]] = np.log10(datadf[settings["lengths_pointer"]])
        settings["lengths_pointer"] = "log_" + settings["lengths_pointer"]
        length_prefix_list.append("Log_")
        logging.info("Using Log10 scaled read lengths.")
        settings["logBool"] = True
    else:
        settings["logBool"] = False
    if settings.get("runtime_until"):
        num_reads_prior = len(datadf)
        datadf = datadf[datadf.start_time < timedelta(hours=settings["runtime_until"])]
        num_reads_post = len(datadf)
        logging.info("Removing {} reads generated after {} hours in the run.".format(
            str(num_reads_prior - num_reads_post),
            str(settings["runtime_until"])))
        settings["filtered"] = True
    if settings.get("downsample"):
        new_size = min(settings["downsample"], len(datadf.index))
        length_prefix_list.append("Downsampled_")
        logging.info("Downsampling the dataset from {} to {} reads".format(
            len(datadf.index), new_size))
        datadf = datadf.sample(new_size)
        settings["filtered"] = True
    if settings.get("percentqual"):
        datadf["quals"] = datadf["quals"].apply(nanomath.phred_to_percent)
        logging.info("Converting quality scores to theoretical percent identities.")
    logging.info("Processed the reads, optionally filtered. {} reads left".format(str(len(datadf))))
    settings["length_prefix"] = ''.join(length_prefix_list)
    return(datadf, settings)
