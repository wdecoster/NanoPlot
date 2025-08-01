#! /usr/bin/env python
# wdecoster

"""
The main purpose of this script is to create plots for long read sequencing data.
Input data can be given as one or multiple of:
-compressed, standard or streamed fastq file
-compressed, standard or streamed fastq file, with
 additional information added by albacore or MinKNOW
-a bam file
-a summary file generated by albacore
"""

from os import path
import logging
import nanoplot.utils as utils
from nanoplot.version import __version__
import sys


def main():
    """
    Organization function
    -setups logging
    -gets inputdata
    -calls plotting function
    """
    settings, args = utils.get_args()
    try:
        import nanoplot.report as report
        from nanoget import get_input
        from nanoplot.filteroptions import filter_and_transform_data
        from nanoplotter.plot import Plot
        utils.make_output_dir(args.outdir)
        import pickle
        utils.init_logs(args)
        # args.format = nanoplotter.check_valid_format(args.format)
        if args.pickle:
            datadf = pickle.load(open(args.pickle, "rb"))
        elif args.feather:
            from nanoget import combine_dfs
            from pandas import read_feather

            datadf = combine_dfs([read_feather(p) for p in args.feather], method="simple").rename(
                columns={"identities": "percentIdentity"}
            )
        else:
            sources = {
                "fastq": args.fastq,
                "bam": args.bam,
                "cram": args.cram,
                "fastq_rich": args.fastq_rich,
                "fastq_minimal": args.fastq_minimal,
                "summary": args.summary,
                "fasta": args.fasta,
                "ubam": args.ubam,
            }
            datadf = get_input(
                source=[n for n, s in sources.items() if s][0],
                files=[f for f in sources.values() if f][0],
                threads=args.threads,
                readtype=args.readtype,
                combine="simple",
                barcoded=args.barcoded,
                huge=args.huge,
                keep_supp=not (args.no_supplementary),
            )
        if args.store:
            pickle.dump(obj=datadf, file=open(settings["path"] + "NanoPlot-data.pickle", "wb"))
        if args.raw:
            datadf.to_csv(
                settings["path"] + "NanoPlot-data.tsv.gz", sep="\t", index=False, compression="gzip"
            )

        settings["statsfile"] = [make_stats(datadf, settings, suffix="", tsv_stats=args.tsv_stats)]
        datadf, settings = filter_and_transform_data(datadf, settings)
        if settings["filtered"]:  # Bool set when filter was applied in filter_and_transform_data()
            settings["statsfile"].append(
                make_stats(
                    datadf[datadf["length_filter"]],
                    settings,
                    suffix="_post_filtering",
                    tsv_stats=args.tsv_stats,
                )
            )

        if args.only_report:
                Plot.only_report = True
                
        if args.barcoded:
            main_path = settings["path"]
            for barc in list(datadf["barcode"].unique()):
                dfbarc = datadf[datadf["barcode"] == barc]
                if len(dfbarc) > 5:
                    logging.info(f"Processing {barc}")
                    settings["title"] = barc
                    settings["path"] = path.join(args.outdir, args.prefix + barc + "_")
                    plots = [report.BarcodeTitle(barc)]
                    plots.extend(make_plots(dfbarc, settings))
                    make_report(plots, settings)
                else:
                    sys.stderr.write(f"Found barcode {barc} less than 5x, ignoring...\n")
                    logging.info(f"Found barcode {barc} less than 5 times, ignoring")
            settings["path"] = main_path
        else:
            plots = make_plots(datadf, settings)
            make_report(plots, settings)
        logging.info("Finished!")
    except Exception as e:
        logging.error(e, exc_info=True)
        print(f"\n\n\nIf you read this then NanoPlot {__version__} has crashed :-(")
        print("Please try updating NanoPlot and see if that helps...\n")
        print("If not, please report this issue at https://github.com/wdecoster/NanoPlot/issues")
        print("If you could include the log file that would be really helpful.")
        print("Thanks!\n\n\n")
        raise


def make_stats(datadf, settings, suffix, tsv_stats=True):
    import nanomath
    statsfile = settings["path"] + "NanoStats" + suffix + ".txt"
    stats_df = nanomath.write_stats(datadfs=[datadf], outputfile=statsfile, as_tsv=tsv_stats)
    logging.info("Calculated statistics")
    if settings["barcoded"]:
        barcodes = list(datadf["barcode"].unique())
        statsfile = settings["path"] + "NanoStats_barcoded" + suffix + ".txt"
        stats_df = nanomath.write_stats(
            datadfs=[datadf[datadf["barcode"] == b] for b in barcodes],
            outputfile=statsfile,
            names=barcodes,
            as_tsv=tsv_stats,
        )
    return stats_df if tsv_stats else statsfile


def make_plots(datadf, settings):
    """
    Call plotting functions from nanoplotter
    settings["lengths_pointer"] is a column in the DataFrame specifying which lengths to use
    """
    import numpy as np
    from scipy import stats
    import nanoplotter
    import nanomath
    color = nanoplotter.check_valid_color(settings["color"])
    colormap = nanoplotter.check_valid_colormap(settings["colormap"])

    plotdict = {type: settings["plots"].count(type) for type in ["kde", "hex", "dot", "pauvre"]}
    if "hex" in settings["plots"]:
        print(
            "WARNING: hex as part of --plots has been deprecated and will be ignored. To get the hex output, rerun with --legacy hex."
        )

    if settings["legacy"]:
        plotdict_legacy = {plot: settings["legacy"].count(plot) for plot in ["kde", "hex", "dot"]}
    else:
        plotdict_legacy = {}
    plots = []

    subdf = utils.subsample_datasets(datadf) if "start_time" in datadf else None
    if settings["N50"]:
        n50 = nanomath.get_N50(np.sort(datadf["lengths"]))
    else:
        n50 = None

    plots.extend(
        nanoplotter.length_plots(
            array=datadf[datadf["length_filter"]]["lengths"].astype("uint64"),
            name="Read length",
            path=settings["path"],
            n50=n50,
            color=color,
            title=settings["title"],
            settings=settings,
        )
    )
    logging.info("Created length plots")
    if "quals" in datadf:
        plots.extend(
            nanoplotter.scatter(
                x=datadf[datadf["length_filter"]][settings["lengths_pointer"].replace("log_", "")],
                y=datadf[datadf["length_filter"]]["quals"],
                legacy=plotdict_legacy,
                names=["Read lengths", "Average read quality"],
                path=settings["path"] + "LengthvsQualityScatterPlot",
                color=color,
                colormap=colormap,
                plots=plotdict,
                title=settings["title"],
                settings=settings,
            )
        )
        if settings["logBool"]:
            plots.extend(
                nanoplotter.scatter(
                    x=datadf[datadf["length_filter"]][settings["lengths_pointer"]],
                    y=datadf[datadf["length_filter"]]["quals"],
                    legacy=plotdict_legacy,
                    names=["Read lengths", "Average read quality"],
                    path=settings["path"] + "LengthvsQualityScatterPlot",
                    color=color,
                    colormap=colormap,
                    plots=plotdict,
                    log=True,
                    title=settings["title"],
                    settings=settings,
                )
            )
        logging.info("Created LengthvsQual plot")
    if "channelIDs" in datadf:
        plots.extend(
            nanoplotter.spatial_heatmap(
                array=datadf["channelIDs"],
                title=settings["title"],
                path=settings["path"] + "ActivityMap_ReadsPerChannel",
                colormap=colormap,
                settings=settings,
            )
        )
        logging.info("Created spatialheatmap for succesfull basecalls.")
    if "start_time" in datadf:
        plots.extend(
            nanoplotter.time_plots(
                df=datadf,
                subsampled_df=subdf,
                path=settings["path"],
                color=color,
                title=settings["title"],
                settings=settings,
            )
        )
        if settings["logBool"]:
            plots.extend(
                nanoplotter.time_plots(
                    df=datadf,
                    subsampled_df=subdf,
                    path=settings["path"],
                    color=color,
                    title=settings["title"],
                    log_length=True,
                    settings=settings,
                )
            )
        logging.info("Created timeplots.")
    if "aligned_lengths" in datadf and "lengths" in datadf:
        plots.extend(
            nanoplotter.scatter(
                x=datadf[datadf["length_filter"]]["aligned_lengths"],
                y=datadf[datadf["length_filter"]]["lengths"],
                legacy=plotdict_legacy,
                names=["Aligned read lengths", "Sequenced read length"],
                path=settings["path"] + "AlignedReadlengthvsSequencedReadLength",
                plots=plotdict,
                color=color,
                colormap=colormap,
                title=settings["title"],
                settings=settings,
            )
        )
        logging.info("Created AlignedLength vs Length plot.")
    if "mapQ" in datadf and "quals" in datadf:
        plots.extend(
            nanoplotter.scatter(
                x=datadf["mapQ"],
                y=datadf["quals"],
                legacy=plotdict_legacy,
                names=["Read mapping quality", "Average basecall quality"],
                path=settings["path"] + "MappingQualityvsAverageBaseQuality",
                color=color,
                colormap=colormap,
                plots=plotdict,
                title=settings["title"],
                settings=settings,
            )
        )
        logging.info("Created MapQvsBaseQ plot.")
        plots.extend(
            nanoplotter.scatter(
                x=datadf[datadf["length_filter"]][settings["lengths_pointer"].replace("log_", "")],
                y=datadf[datadf["length_filter"]]["mapQ"],
                legacy=plotdict_legacy,
                names=["Read length", "Read mapping quality"],
                path=settings["path"] + "MappingQualityvsReadLength",
                color=color,
                colormap=colormap,
                plots=plotdict,
                title=settings["title"],
                settings=settings,
            )
        )
        if settings["logBool"]:
            plots.extend(
                nanoplotter.scatter(
                    x=datadf[datadf["length_filter"]][settings["lengths_pointer"]],
                    y=datadf[datadf["length_filter"]]["mapQ"],
                    legacy=plotdict_legacy,
                    names=["Read length", "Read mapping quality"],
                    path=settings["path"] + "MappingQualityvsReadLength",
                    color=color,
                    colormap=colormap,
                    plots=plotdict,
                    log=True,
                    title=settings["title"],
                    settings=settings,
                )
            )
        logging.info("Created Mapping quality vs read length plot.")
    if "percentIdentity" in datadf:
        minPID = np.percentile(datadf["percentIdentity"], 1)
        if "aligned_quals" in datadf:
            plots.extend(
                nanoplotter.scatter(
                    x=datadf["percentIdentity"],
                    y=datadf["aligned_quals"],
                    legacy=plotdict_legacy,
                    names=["Percent identity", "Average Base Quality"],
                    path=settings["path"] + "PercentIdentityvsAverageBaseQuality",
                    color=color,
                    colormap=colormap,
                    plots=plotdict,
                    stat=stats.pearsonr if not settings["hide_stats"] else None,
                    minvalx=minPID,
                    title=settings["title"],
                    settings=settings,
                )
            )
            logging.info("Created Percent ID vs Base quality plot.")
        plots.extend(
            nanoplotter.scatter(
                x=datadf[datadf["length_filter"]][settings["lengths_pointer"].replace("log_", "")],
                y=datadf[datadf["length_filter"]]["percentIdentity"],
                legacy=plotdict_legacy,
                names=["Aligned read length", "Percent identity"],
                path=settings["path"] + "PercentIdentityvsAlignedReadLength",
                color=color,
                colormap=colormap,
                plots=plotdict,
                stat=stats.pearsonr if not settings["hide_stats"] else None,
                minvaly=minPID,
                title=settings["title"],
                settings=settings,
            )
        )
        if settings["logBool"]:
            plots.extend(
                nanoplotter.scatter(
                    x=datadf[datadf["length_filter"]][settings["lengths_pointer"]],
                    y=datadf[datadf["length_filter"]]["percentIdentity"],
                    legacy=plotdict_legacy,
                    names=["Aligned read length", "Percent identity"],
                    path=settings["path"] + "PercentIdentityvsAlignedReadLength",
                    color=color,
                    colormap=colormap,
                    plots=plotdict,
                    stat=stats.pearsonr if not settings["hide_stats"] else None,
                    log=True,
                    minvaly=minPID,
                    title=settings["title"],
                    settings=settings,
                )
            )

        plots.append(
            nanoplotter.dynamic_histogram(
                array=datadf["percentIdentity"],
                name="percent identity",
                path=settings["path"] + "PercentIdentityHistogram",
                title=settings["title"],
                color=color,
                settings=settings,
            )
        )
        logging.info("Created Percent ID vs Length plot")
    return plots


def make_report(plots, settings):
    """
    Creates a fat html report based on the previously created files
    plots is a list of Plot objects defined by a path and title
    statsfile is the file to which the stats have been saved,
    which is parsed to a table (rather dodgy) or nicely if it's a pandas/tsv
    """
    logging.info("Writing html report.")
    from nanoplot import report

    html_content = [
        '<body class="grid">',
        report.html_toc(plots, filtered=settings["filtered"]),
        report.html_stats(settings),
        report.html_plots(plots),
        report.run_info(settings) if settings["info_in_report"] else "",
        "</main></body></html>",
    ]
    with open(settings["path"] + "NanoPlot-report.html", "w") as html_file:
        html_file.write(report.html_head + "\n".join(html_content))


if __name__ == "__main__":
    main()
