import os
import sys
import subprocess
import pathlib
import collections
from itertools import groupby
from Bio import SeqIO
import matplotlib.pyplot as plt

__author__ = 'Colin Anthony'


def try_except_exit_on_fail(cmd):
    try:
        subprocess.call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(e)
        sys.exit("exiting")


def try_except_continue_on_fail(cmd):
    try:
        subprocess.call(cmd, shell=True)
        return True
    except subprocess.CalledProcessError as e:
        print(e)
        return False

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
        new_key = k.replace(" ", "_") + "_" + str(i).zfill(6)
        dct[new_key] = v.upper()

    return dct


def gather_fastqs(fastq_path, run_name, max_len, min_len):

    fastq_outpath = fastq_path.parent
    fastq_name = f"{run_name}_all.fastq"
    output_fastq = pathlib.Path(fastq_outpath, fastq_name)

    all_fastqs = list(fastq_path.glob("*.fastq"))

    with open(output_fastq, 'w') as handle:
        for fastq in all_fastqs:
            try:
                for record in SeqIO.parse(open(fastq), "fastq"):
                    seq_len = len(record.seq)
                    if seq_len > max_len or seq_len < min_len:
                        continue
                    else:
                        SeqIO.write([record], handle, "fastq")
            except ValueError as e:
                print("Failed on fastq file:", fastq, "\n", e, "\n", "Continuing with next fastq file")
                continue
    if output_fastq.is_file():
        return True
    else:
        return False


def filter_length(fastq, outfile, max_len, min_len):
    """
    filter fastq by length
    :param fastq: (str) the path and name of the fastq file
    :param outfile: (str) the path and name of the output fastq file
    :param max_len: (int) max sequence length
    :param min_len: (int) min sequence length
    :return: (bool) True if outfile is not empty else False
    """

    not_empty = False
    with open(outfile, 'w') as handle:
        try:
            for record in SeqIO.parse(open(fastq), "fastq"):
                seq_len = len(record.seq)
                if seq_len > max_len or seq_len < min_len:
                    continue
                else:
                    not_empty = True
                    SeqIO.write([record], handle, "fastq")
        except ValueError as e:
            print("Failed on fastq file:", fastq, "\n", e, "\n", "Continuing with next fastq file")
            return None

    return not_empty

def filter_length_trim_seq(fastq, outfile, max_len, min_len, left_trim, right_trim):
    """
    filter fastq by length and trim nt from left and right
    :param fastq: (str) the path and name of the fastq file
    :param outfile: (str) the path and name of the output fastq file
    :param max_len: (int) max sequence length
    :param min_len: (int) min sequence length
    :param left_trim: (int) number of nt to trim from left
    :param right_trim: (int) number of nt to trim from right
    :return: (bool) True if outfile is not empty else False
    """
    left_trim_forrmat = left_trim + 1
    right_trim_format = -right_trim
    not_empty = False
    with open(outfile, 'w') as handle:
        try:
            for record in SeqIO.parse(open(fastq), "fastq"):
                seq_len = len(record.seq)
                if seq_len > max_len or seq_len < min_len:
                    continue
                else:
                    not_empty = True
                    trimmed_record = record[left_trim_forrmat:right_trim_format]
                    SeqIO.write([trimmed_record], handle, "fastq")
        except ValueError as e:
            print("Failed on fastq file:", fastq, "\n", e, "\n", "Continuing with next fastq file")
            return None

    return not_empty

def d_freq_lists_pos(dna_list, n, positional_depth):
    """
    calculate base frequencies from list of sequences (aligned) using positional depth
    :param dna_list: (list) a list of DNA sequences
    :param n: (int) length of alignment (ie number of positions)
    :param positional_depth: (dict) key = str(position), value = int(depth)
    :return: (dict) a dictionary of the frequency for each base, for each site in the alignment
    """
    counts_dict = {'A': [0]*n, 'C': [0]*n, 'G': [0]*n, 'T': [0]*n, '-': [0]*n, "N": [0]*n}
    bases = ["A", "C", "G", "T", "-"]
    depth_dict = {'non_gap': [0] * n, 'gap': [0] * n}
    for seq in dna_list:
        for index, dna in enumerate(seq):
            if dna.upper() not in bases:
                counts_dict["N"][index] += 1
            else:
                counts_dict[dna.upper()][index] += 1
            if dna == "-":
                depth_dict["gap"][index] += 1
            else:
                depth_dict["non_gap"][index] += 1

    total_seqs = len(dna_list)
    for base, countslist in counts_dict.items():
        for position, cnt in enumerate(countslist):
            position_lookup = str(position +1)
            try:
                position_depth = positional_depth[position_lookup]
            except IndexError:
                print(f"Positional depth could not be calculated from primer-pair depth for position {position_lookup}")
                if base != '-':
                    position_depth = depth_dict["non_gap"][position]
                else:
                    position_depth = 0
            if position_depth == 0:
                frq = 0
            else:
                # adjust count for gaps for only the primer pair coverage
                if base == "-":
                    cnt = position_depth - (total_seqs - cnt)
                frq = round((cnt/position_depth*100), 4)
            countslist[position] = frq
        counts_dict[base] = countslist

    return counts_dict, depth_dict


def consensus_maker(d, positional_depth, min_depth, use_gaps):
    """
    Create a consensus sequence from an alignment
    :param d: (dict) dictionary of an alignment (key = seq name (str): value = aligned sequence (str))
    :param positional_depth: (dict) key = str(position), value = int(depth)
    :param min_depth: (int) the minimum depth required to call a base in the consensus (otherwise called as "!"
    :param use_gaps (bool) use gap characters when making consensus
    :return: (str) the consensus sequence
    """
    seq_list = []
    for names, seq in d.items():
        seq_list.append(seq)
    if not seq_list:
        raise IndexError

    seq_length = len(seq_list[0])
    master_profile, depth_profile = d_freq_lists_pos(seq_list, seq_length, positional_depth)
    consensus = ""
    degen = {('A', 'G'): 'R', ('C', 'T'): 'Y', ('A', 'C'): 'M', ('G', 'T'): 'K', ('C', 'G'): 'S', ('A', 'T'): 'W',
             ('A', 'C', 'T'): 'H', ('C', 'G', 'T'): 'B', ('A', 'C', 'G'): 'V', ('A', 'G', 'T'): 'D',
             ('A', 'C', 'G', 'T'): 'N'}

    for position in range(seq_length):
        position_lookup = str(position+1)
        if positional_depth[position_lookup] <= min_depth:
            consensus += str("N")
        else:
            if use_gaps:
                dct = {base: master_profile[base][position] for base in ['A', 'C', 'G', 'T', 'N', "-"]}
            else:
                dct = {base: master_profile[base][position] for base in ['A', 'C', 'G', 'T', 'N']}
            # get the base with the highest frequency value
            base_with_max_freq = max(dct, key=dct.get)
            # get the highest frequency value
            max_freq = dct[base_with_max_freq]
            # if multiple bases share the max frequency make a list of them for degeneracy code lookup
            if use_gaps:
                most_freq_bases = list(sorted(base for base in ['A', 'C', 'G', 'T', "-"] if dct[base] == max_freq))
            else:
                most_freq_bases = list(sorted(base for base in ['A', 'C', 'G', 'T'] if dct[base] == max_freq))

            if len(most_freq_bases) == 1:
                consensus += str(base_with_max_freq)
            else:
                if '-' in most_freq_bases:
                    most_freq_bases.remove("-")
                if len(most_freq_bases) == 1:
                    consensus += str(base_with_max_freq)
                else:
                    most_freq_bases = tuple(most_freq_bases)
                    degen_char = str(degen[most_freq_bases])
                    consensus += degen_char

    return consensus, depth_profile


def autolabel(barchart, axis, primer_depth):
    """
    Attach a text label above each bar displaying its height
    """
    for i, bar in enumerate(barchart):
        height = bar.get_height()
        label = primer_depth[i]

        axis.text(bar.get_x() + bar.get_width()/2., 1.05*height,
                f'{label}', ha='center', va='bottom')


def plot_primer_depth(primer_pairs, primer_depth, percent_primers_depth, sample_name, outfile):

    y_vals = percent_primers_depth

    fig, ax = plt.subplots()
    ax.set_ylabel('Sequencing depth \n(% of max primer pair depth)')
    ax.set_xlabel('Primer pair')
    ax.set_title(sample_name)

    ax.set_ylim(bottom=0, top=121)

    # plot bar graph
    bar_plot = plt.bar(primer_pairs, y_vals)
    plt.xticks(rotation='vertical')

    # add value to each bar
    autolabel(bar_plot, ax, primer_depth)

    w = 6.8
    h = 4
    f = plt.gcf()
    f.set_size_inches(w, h)

    plt.savefig(outfile, ext="png", dpi=300, facecolor="white", bbox_inches="tight")
    plt.close()


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
    plt.close()


def rename_fasta(fasta_file_name_path, sample_name, cons_type):
    fasta_d = fasta_to_dct(fasta_file_name_path)
    os.unlink(fasta_file_name_path)
    with open(fasta_file_name_path, 'w') as fh:
        for seq_name, seq in fasta_d.items():
            new_name = f"{sample_name}_{cons_type}"
            fh.write(f">{new_name}\n{seq}\n")


def cat_sample_names(barcode, run_name):
    if barcode != '':
        file_name = f"{run_name}_{barcode}.fastq"
    else:
        file_name = " "

    return file_name