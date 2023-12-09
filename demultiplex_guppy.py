import argparse
import pathlib
from src.misc_functions import try_except_continue_on_fail


__author__ = 'Philippe Selhorst'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass

def main(inpath, guppy_path, outpath, one_end):
    # force absolute file paths
    both_ends=''
    if not one_end:
        both_ends = '--require_barcodes_both_ends'
        print(f'{both_ends}\n')
    if one_end:
        print('--require_barcodes_one_end\n')
    inpath = pathlib.Path(inpath).absolute()
    outpath = pathlib.Path(outpath).absolute()
    guppy_path = pathlib.Path(guppy_path).absolute()
    guppy_demultiplexer = pathlib.Path(guppy_path, "guppy_barcoder")
    gpu_settings = f"-x 'auto'"

    guppy_demux_cmd = f"{str(guppy_demultiplexer)} -i {inpath} -s {outpath} " \
                      f"--enable_trim_barcodes -t 24 --num_barcoding_threads 16 " \
                      f"--records_per_fastq 0 {gpu_settings} {both_ends}"

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
    parser.add_argument('-oe', '--one_end', type=str, default=None, required=False,
                        help='The path for the outfile')

    args = parser.parse_args()
    inpath = args.inpath
    guppy_path = args.guppy_path
    outpath = args.outpath
    one_end = args.one_end

    main(inpath, guppy_path, outpath, one_end)
