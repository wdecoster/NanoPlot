import pandas as pd
from argparse import ArgumentParser


def main():
    parser = ArgumentParser(
        description="Add barcode_arrangement to sequencing_summary.txt using barcoding_summary.txt")
    parser.add_argument("sequencing_summary",
                        help="sequencing_summary.txt from guppy, can be compressed.")
    parser.add_argument("barcoding_summary",
                        help="barcoding_summary.txt from guppy, can be compressed.")
    args = parser.parse_args()
    bc = pd.read_csv(args.barcoding_summary, sep="\t", usecols=['read_id', 'barcode_arrangement'])
    df = pd.read_csv(args.sequencing_summary, sep="\t")
    df.join(bc.set_index('read_id'), on='read_id') \
      .to_csv(path_or_buf="barcoded_sequencing_summary.txt.gz",
              sep="\t",
              compression="gzip",
              index=False)


if __name__ == '__main__':
    main()
