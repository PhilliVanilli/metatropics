import collections
from itertools import groupby
import pathlib
import argparse
import sys


__author__ = 'Colin Anthony'


def py3_fasta_iter(fasta_name):
    """
    modified from Brent Pedersen: https://www.biostars.org/p/710/#1412
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(str(fasta_name), 'r')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header_str = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (header_str, seq)


def fasta_to_dct(file_name):
    """
    :param file_name: The fasta formatted file to read from.
    :return: a dictionary of the contents of the file name given. Dictionary in the format:
    {sequence_id: sequence_string, id_2: sequence_2, etc.}
    """
    dct = collections.defaultdict(str)
    my_gen = py3_fasta_iter(file_name)
    for i, (k, v) in enumerate(my_gen):
        # resolve for duplicate sequence headers
        new_key = k.replace(" ", "_")
        dct[new_key] = v.upper()

    return dct


def main(infile, outpath, ref_name, scheme_name):

    seqs_d = fasta_to_dct(infile)

    try:
        ref_sequence = seqs_d[ref_name]
    except KeyError as e:
        print("Reference sequence name not found in fasta file\nExiting\n")
        sys.exit(e)

    del seqs_d[ref_name]

    ref_bed = pathlib.Path(outpath, f"{scheme_name}.scheme.bed")
    ref_fasta = pathlib.Path(outpath, f"{scheme_name}.reference.fasta")

    with open(ref_fasta, 'w') as fh:
        fh.write(f">{ref_name}\n{ref_sequence}\n")

    collected_d = collections.defaultdict(list)
    for seq_name, seq in seqs_d.items():
        first = True
        primer_start = None
        primer_end = None
        for indx, pos in enumerate(seq):
            if pos != "-" and first:
                primer_start = indx
                first = False

            if pos != "-" and seq[indx + 1] == "-":
                primer_end = indx
        if primer_start is None or primer_end is None:
            print(f"could not find primer entry {seq_name}\n skipping it")
            continue
        collected_d[seq_name] = [primer_start, primer_end, seq]

    with open(ref_bed, 'w') as fh1:
        fh1.write("genome\tstart\tend\tPrimer_ID\tnumber\n")

        for primer, (primer_start, primer_end, seq) in collected_d.items():
            num = primer.split("_")[1]
            len_primer = len(seq.replace("-", ""))
            ref_start = len(ref_sequence[:primer_start].replace("-", ""))
            ref_end = ref_start + len_primer

            fh1.write(f"{ref_name}\t{ref_start}\t{ref_end}\t{primer}\t{num}\n")

    print("done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="make bed file from primer to ref alignment",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-in", "--infile", default=argparse.SUPPRESS, type=str,
                        help="The path and name of the fasta file with reference and primer sequences aligned",
                        required=True)
    parser.add_argument("-o", "--outpath", default=argparse.SUPPRESS, type=str,
                        help="The path to where the output bed files will be written", required=True)
    parser.add_argument("-r", "--ref_name", default=argparse.SUPPRESS, type=str,
                        help="The name of the reference sequence in the fasta file: ", required=True)
    parser.add_argument("-s", "--scheme_name", default=argparse.SUPPRESS, type=str,
                        help="The name of the reference scheme. id: ChikAsian400", required=True)

    args = parser.parse_args()

    infile = args.infile
    outpath = args.outpath
    ref_name = args.ref_name
    scheme_name = args.scheme_name

    main(infile, outpath, ref_name, scheme_name)
