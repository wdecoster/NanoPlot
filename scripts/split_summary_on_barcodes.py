import pandas as pd
from argparse import ArgumentParser


def main():
    parser = ArgumentParser(description="Split sequencing_summary.txt using barcoding_summary.txt")
    parser.add_argument("sequencing_summary",
                        help="sequencing_summary.txt from guppy, can be compressed.")
    parser.add_argument("barcoding_summary",
                        help="barcoding_summary.txt from guppy, can be compressed.")
    args = parser.parse_args()
    bc = pd.read_csv(args.barcoding_summary, sep="\t", usecols=['read_id', 'barcode_arrangement'])
    df = pd.read_csv(args.sequencing_summary, sep="\t")
    for barc in bc["barcode_arrangement"].unique():
        df[df["read_id"].isin(bc.loc[bc["barcode_arrangement"] == barc, 'read_id'])] \
            .to_csv("sequencing_summary_{}.txt".format(barc), sep="\t", index=False)


if __name__ == '__main__':
    main()
