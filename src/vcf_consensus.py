import argparse
import pathlib
import vcfpy
import vcf
import subprocess
import os
from collections import OrderedDict


__author__ = 'Colin Anthony'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def main(infile, min_depth, outpath):
    # force absolute file paths
    infile = pathlib.Path(infile).absolute()
    outpath = pathlib.Path(outpath).absolute()
    cons_outfile = pathlib.Path(outpath, infile.stem + "_vcf_consensus.fasta")
    depth_qual_outfile = pathlib.Path(outpath, infile.stem + "_depth_qual.csv")
    headers = "position,seq_depth,quality_score\n"
    sample_name = infile.stem.replace("_bcftools", "")

    cmd = f"bcftools view -Ov -o {infile}.decomplressed.vcf {infile}"
    subprocess.call(cmd, shell=True)
    ref_seq = OrderedDict()
    cons_seq = OrderedDict()
    first = True
    csv_string = ""
    ref_name = ""
    real_del_list = []
    for record in vcfpy.Reader.from_path(str(infile) + ".decomplressed.vcf"):
        pos = record.POS
        if first:
            ref_name = record.CHROM
            cons_seq["0"] = "N" * (pos - 1)
            ref_seq["0"] = "N" * (pos - 1)
            first = False

        ref_base = record.REF
        ref_seq[pos] = ref_base
        alt_subs = record.ALT
        qual =record.QUAL
        inf = record.INFO
        depth = inf["DP"]

        if not record.is_snv():
            indel_freq = record.INFO["IMF"]
            if indel_freq >= 0.6:
                del_len = len(ref_base)
                real_del_list.extend(list(range(pos, pos + del_len, 1)))
            continue

        csv_string += f"{pos},{depth},{qual}\n"

        if len(alt_subs) >= 1 and depth >= min_depth and qual >= 30:
            alt_base = record.ALT[0].value
        else:
            alt_base = ref_base
        cons_seq[pos] = alt_base
        # cons_seq += alt_base
    consensus = ""
    reference = ""

    for position, base in cons_seq.items():
        if position in real_del_list:
            continue
        consensus += base.upper()
    for position, base in ref_seq.items():
        reference += base.upper()
    # with cons_outfile.open("w") as fh:
    #     fh.write(f">{sample_name}\n{consensus}\n")

    with depth_qual_outfile.open('w') as fh:
        fh.write(headers)
        fh.write(csv_string)

    # with open("ref.fasta", "w") as fh:
    #     fh.write(f">{ref_name}\n{reference}\n")

    print("done")
    # os.unlink(str(infile) + ".decomplressed.vcf")

    return depth_qual_outfile


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=Formatter)

    parser.add_argument('-in', '--infile', type=str, default=None, required=True,
                        help='The path and name of the infile')
    parser.add_argument('-min', '--mind_depth', type=str, default=100, required=False,
                        help='The path and name of the infile')
    parser.add_argument('-o', '--outpath', type=str, default=None, required=True,
                        help='The path for the outfile')

    args = parser.parse_args()
    infile = args.infile
    mind_depth = args.mind_depth
    outpath = args.outpath

    main(infile, mind_depth, outpath)
