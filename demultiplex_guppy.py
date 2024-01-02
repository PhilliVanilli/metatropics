import argparse
import pathlib
from src.misc_functions import try_except_continue_on_fail


__author__ = 'Philippe Selhorst'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass

def main(inpath, guppy_path, outpath, barcodes, one_end):
    # force absolute file paths
    ends = ""
    if not one_end:
        ends = "--require_barcodes_both_ends"
        print(f"{ends}")
    if one_end:
        print("--require_barcodes_one_end\n")
    inpath = pathlib.Path(inpath).absolute()
    outpath = pathlib.Path(outpath).absolute()
    guppy_path = pathlib.Path(guppy_path).absolute()
    guppy_demultiplexer = pathlib.Path(guppy_path, "guppy_barcoder")
    gpu_settings = f"-x 'auto'"
    barcode_args = f"--barcode_kits {barcodes}"

    guppy_demux_cmd = f"{str(guppy_demultiplexer)} -i {inpath} -s {outpath} " \
                      f"--enable_trim_barcodes -t 48 --num_barcoding_threads 16 " \
                      f"--records_per_fastq 0 {gpu_settings} {ends} {barcode_args}"

    run = try_except_continue_on_fail(guppy_demux_cmd)

    if run:
        print("Demultiplexing completed\n")
    else:
        print("Demultiplexing failed")

    return run


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to wrap guppy demultiplex commands',
                                     formatter_class=Formatter)

    parser.add_argument('-in', '--inpath', type=str, default=None, required=True,
                        help='The path to the fastq folder')
    parser.add_argument('-p', '--guppy_path', type=str, default=None, required=True,
                        help='The path to guppy executable')
    parser.add_argument('-o', '--outpath', type=str, default=None, required=True,
                        help='The path for the outfile')
    parser.add_argument("-bc", "--barcodes", type=str, choices=["CUST","SQK-NBD114-24"], required=True,
                        help="Specify barcodes used for demultiplexing")
    parser.add_argument("-oe", "--one_end", default=False, type=str, action="store_true",
                        help="use reads if they have barcode on only one end, this increases the amount of data yet increases probability of misclassification")
    args = parser.parse_args()
    inpath = args.inpath
    guppy_path = args.guppy_path
    outpath = args.outpath
    barcodes = args.barcodes
    one_end = args.one_end

    main(inpath, guppy_path, outpath, barcodes, one_end)
