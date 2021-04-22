#!/usr/bin/env python
from Bio import SeqIO
import sys
import vcf
import subprocess
from collections import defaultdict
import os.path
import argparse
import pathlib
import matplotlib.pyplot as plt


def collect_depths(bamfile):
    if not os.path.exists(bamfile):
        raise SystemExit("bamfile %s doesn't exist" % (bamfile,))

    print((bamfile, sys.stderr))

    p = subprocess.Popen(['samtools', 'depth', bamfile], stdout=subprocess.PIPE)
    out, err = p.communicate()
    depths = defaultdict(dict)
    # print(out[0])
    # input("dfdf")
    for ln in out.decode('utf-8').split("\n"):
        if ln:
            contig, pos, depth = ln.split("\t")
            depths[contig][int(pos)] = int(depth)

    return depths


def plot_depth(depth_list, sample_name, outfile):

    x_vals = [x for x in range(len(depth_list))]
    fig, ax = plt.subplots()
    ax.set_ylabel('Sequencing depth')
    ax.set_xlabel('Sequence position')
    ax.set_title(sample_name)

    plt.plot(x_vals, depth_list)

    w = 6.8
    h = 4
    f = plt.gcf()
    f.set_size_inches(w, h)

    plt.savefig(outfile, ext="png", dpi=300, facecolor="white", bbox_inches="tight")


def plot_qual(qual_list, sample_name, outfile):

    x_vals = [x for x in range(len(qual_list))]
    fig, ax = plt.subplots()
    ax.set_ylabel('Sequencing quality')
    ax.set_xlabel('Sequence position')
    ax.set_title(sample_name)

    plt.plot(x_vals, qual_list)

    plt.savefig(outfile, ext="png", dpi=300, facecolor="white", bbox_inches="tight")


def main(reference, vcf_file, bamfile, sample_name, outpath):

    if bamfile:
        bamfile = pathlib.Path(bamfile).absolute()
        reference = pathlib.Path(reference).absolute()
        outfile_depth = pathlib.Path(outpath, sample_name + "_sequencing_depth.png")

        depths = collect_depths(bamfile)
        seq = list(SeqIO.parse(open(str(reference)), "fasta"))[0]
        cons = list(seq.seq)

        seq_depth_by_pos = []
        for n, c in enumerate(cons):
            try:
                depth = depths[seq.id][n+1]

            except KeyError:
                depth = 0
            seq_depth_by_pos.append(depth)

        plot_depth(seq_depth_by_pos, sample_name, outfile_depth)

    if vcf_file:
        vcffile = pathlib.Path(vcf_file).absolute()
        vcf_reader = vcf.Reader(open(str(vcffile), 'r'))
        outfile_qual = pathlib.Path(vcffile.parent, sample_name + "_sequencing_qual.png")
        seq_qual_by_pos = []
        for record in vcf_reader:
            seq_qual_by_pos.append(record.QUAL)
        plot_qual(seq_qual_by_pos, sample_name, outfile_qual)

    if not bamfile and not vcf_file:
        print("no files given, no plots made")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="make consensus from vcf and bam file",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-r", "--reference", type=str, default=argparse.SUPPRESS,
                        help="The reference genome and primer scheme to use", required=True)
    parser.add_argument("-o", "--outpath", type=str, default=False,
                        help="The path where the output will be written", required=False)
    parser.add_argument("-v", "--vcf_file", type=str, default=False,
                        help="The path and name of the vcf file", required=False)
    parser.add_argument("-b", "--bam_file", default=False, type=str,
                        help="The path and name of the sorted, trimmed bam file", required=False)
    parser.add_argument("-n", "--sample_name", type=str, default=argparse.SUPPRESS,
                        help="The sample name", required=True)

    args = parser.parse_args()
    vcf_file = args.vcf_file
    bam_file = args.bam_file
    reference = args.reference
    sample_name = args.sample_name
    outpath = args.outpath

    main(reference, vcf_file, bam_file, sample_name, outpath)
