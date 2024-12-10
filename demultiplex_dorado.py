import argparse
import pathlib
import os
import sys
from src.misc_functions import try_except_continue_on_fail

__author__ = 'Philippe Selhorst'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass

def main(inpath, dorado_path, outpath, barcodes, one_end):
    # force absolute file paths
    ends = ""
    if not one_end:
        ends = "--barcode-both-ends"
        print(f"{ends}\n")
    if one_end:
        print("--barcode_one_end\n")
    inpath = pathlib.Path(inpath).absolute()
    pre_existing_files = list(inpath.glob("*fastq*"))
    if pre_existing_files:
        answer = input("Previous demultiplexed files exist, overwrite (y/n)?")
        if answer == 'n':
            sys.exit("Keeping demultiplexed files")
        else:
            for file in pre_existing_files:
                os.remove(file)
    # calls_file = pathlib.Path(inpath,"calls.fastq")

    outpath = pathlib.Path(outpath).absolute()
    dorado_path = pathlib.Path(dorado_path).absolute()
    dorado_demultiplexer = pathlib.Path(dorado_path, "dorado demux")
    # gpu_settings = f"-x 'auto'"

    if barcodes == 'CUST':
        custom_arr = pathlib.Path(dorado_path, "barcode_arrs_cust_dorado.toml")
        custom_seq = pathlib.Path(dorado_path, "barcodes_cust.fastq")
        dorado_demux_cmd = f"{str(dorado_demultiplexer)} -o {outpath} " \
                           f"--emit-fastq {ends} --barcode-arrangement {custom_arr} --barcode-sequences {custom_seq} {inpath}"

    else:
        dorado_demux_cmd = f"{str(dorado_demultiplexer)} --kit-name {barcodes} -o {outpath} " \
                           f"--emit-fastq {ends} {inpath}"

    run = try_except_continue_on_fail(dorado_demux_cmd)

    if run:
        print("Demultiplexing completed\n")
    else:
        print("Demultiplexing failed")

    return run


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to wrap dorado demultiplex commands',
                                     formatter_class=Formatter)

    parser.add_argument('-in', '--inpath', type=str, default=None, required=True,
                        help='The path to the fastq folder')
    parser.add_argument('-p', '--dorado_path', type=str, default=None, required=True,
                        help='The path to dorado executable')
    parser.add_argument('-o', '--outpath', type=str, default=None, required=True,
                        help='The path for the outfile')
    parser.add_argument("-bc", "--barcodes", type=str, choices=["CUST","SQK-NBD114-24"], required=True,
                        help="Specify barcodes used for demultiplexing")
    parser.add_argument("-oe", "--one_end", default=False, type=str, action="store_true",
                        help="use reads if they have barcode on only one end, this increases the amount of data yet increases probability of misclassification")
    args = parser.parse_args()
    inpath = args.inpath
    dorado_path = args.dorado_path
    outpath = args.outpath
    barcodes = args.barcodes
    one_end = args.one_end

    main(inpath, dorado_path, outpath, barcodes, one_end)
