import os
import argparse
import pathlib
from pycoQC.Fast5_to_seq_summary import Fast5_to_seq_summary
import pandas as pd


__author__ = 'Colin Anthony'


def main(fast5_dir, outpath):
    fast5_dir = pathlib.Path(fast5_dir).absolute()
    outpath = pathlib.Path(outpath).absolute()

    print("Creating sequencing_summary.txt file for each fast5 folder\n")
    individual_dir = pathlib.Path(fast5_dir).glob("*")
    fast5_individual_dirs = [f for f in individual_dir if f.is_dir()]
    for folder in fast5_individual_dirs:
        file_num = str(folder.parts[-1]).zfill(3)
        outfile = pathlib.Path(outpath, f"sequencing_summary_{file_num}.txt")
        Fast5_to_seq_summary(fast5_dir=folder, seq_summary_fn=outfile, threads=6, verbose_level=1, include_path=True)

    print("Creating filename column for sequencing_summary.txt file\n")
    for summary_file in pathlib.Path(outpath).glob("sequencing_summary_*.txt"):

        file_df = pd.read_csv(summary_file, sep=None, engine="python")
        file_df["filename"] = file_df["path"].apply(lambda x: x.split("/")[-1])
        os.unlink(str(summary_file))
        file_df.to_csv(summary_file, sep="\t", index=False)

    print("sequencing_summary.txt file creation complete")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="generate sequence_summary.txt file from fast5 file",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--fast5_dir", default=argparse.SUPPRESS, type=str,
                        help="The path to the 'fast5' directory ", required=True)
    parser.add_argument("-o", "--outpath", default=argparse.SUPPRESS, type=str,
                        help="The path to where the sequence_summary file will be written", required=True)
    args = parser.parse_args()

    fast5_dir = args.fast5_dir
    outpath = args.outpath

    main(fast5_dir, outpath)
