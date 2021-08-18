import argparse
import pathlib
import os
import shutil
import csv
from statistics import mean

from src.misc_functions import try_except_continue_on_fail
from src.misc_functions import try_except_exit_on_fail
from src.misc_functions import consensus_maker
from src.misc_functions import fasta_to_dct
from src.misc_functions import plot_depth

__author__ = 'Colin Anthony'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def main(infile, log_file, chosen_ref_file, threads,
         min_depth, use_gaps):

    # force absolute file paths
    sample_fastq = pathlib.Path(infile).absolute()
    script_dir = pathlib.Path(__file__).absolute().parent
    if not sample_fastq.is_file():
        print(f"\nCould not find the concatenated sample fastq file: {sample_fastq}\nskipping sample")
        with open(log_file, "a") as handle:
            handle.write(f"\nCould not find the concatenated sample fastq file: {sample_fastq}\nskipping sample")
        return False

    # set input and output file paths
    sample_name = pathlib.Path(sample_fastq).stem
    sample_dir = pathlib.Path(sample_fastq).parent
    project_dir = sample_dir.parent.parent
    seq_folder = pathlib.Path(project_dir, "seq_files")
    seq_folder.mkdir(mode=0o777, parents=True, exist_ok=True)
    sample_viruses_file = pathlib.Path(seq_folder, sample_name + "_viruses.fasta")
    plot_folder = pathlib.Path(project_dir, "seq_depth_plots")

    # initialize the fasta files with viruses per sample
    reference_d = fasta_to_dct(chosen_ref_file)
    for k,v in reference_d.items():
        ref_name = k[0:-7]
        ref_seq = v.replace('-', '').lower()
        with open(sample_viruses_file,'a') as handle:
            handle.write(f">{ref_name}\n{ref_seq}\n")
            handle.close()

    # iterate mapping over the references
    for k, v in reference_d.items():

        # generate virus dir and individual reference files
        ref_name = k[0:-7]
        print(ref_name)
        ref_seq = v.replace('-', '')
        reference_slice = f"{ref_name}:1-{len(ref_seq)}"
        single_ref_file = pathlib.Path(sample_dir, ref_name + ".fasta")
        with open(single_ref_file, 'w') as fh:
            fh.write(f">{ref_name}\n{ref_seq}\n")
        virus_dir = pathlib.Path(sample_dir, ref_name)
        virus_dir.mkdir(mode=0o777, parents=True, exist_ok=True)
        os.chdir(virus_dir)

        # set input and output file paths
        sam_file = pathlib.Path(virus_dir, sample_name + ".sam")
        bam_file = pathlib.Path(virus_dir, sample_name + "_mapped.bam")
        sorted_bam_file = pathlib.Path(virus_dir, sample_name + "_sorted.bam")
        depth_file = pathlib.Path(virus_dir, sample_name + "_depth.tsv")
        msa_fasta = pathlib.Path(virus_dir, sample_name + "_msa_from_bam_file.fasta")
        msa_cons = pathlib.Path(virus_dir, sample_name + "_msa_consensus.fasta")


        # # run read mapping using minimap
        # print(f"\nRunning: minimap2 read mapping")
        # minimap2_cmd = f"minimap2 -a -Y -t 8 -x ava-ont {single_ref_file} {sample_fastq} -o {sam_file} " \
        #                f"2>&1 | tee -a {log_file}"
        # print("\n", minimap2_cmd, "\n")
        # with open(log_file, "a") as handle:
        #     handle.write(f"\nRunning: minimap read mapping\n")
        #     handle.write(f"{minimap2_cmd}\n")
        # run = try_except_continue_on_fail(minimap2_cmd)
        # if not run:
        #     return False

        # run read mapping using bwa
        make_index_cmd = f"bwa index {single_ref_file}"
        with open(log_file, "a") as handle:
            handle.write(f"\n{make_index_cmd}\n")

        try_except_exit_on_fail(make_index_cmd)

        print(f"\nrunning: bwa read mapping\n")
        bwa_cmd = f"bwa mem -t {threads} -x ont2d {single_ref_file} {sample_fastq} -o {sam_file} " \
                  f"2>&1 | tee -a {log_file}"
        print("\n", bwa_cmd, "\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nrunning: bwa read mapping\n")
            handle.write(f"{bwa_cmd}\n")
        run = try_except_continue_on_fail(bwa_cmd)
        if not run:
            return False


        # convert sam to bam
        print(f"\nRunning: sam to bam conversion of mapped file")
        sam_bam_cmd = f"samtools view -bS {sam_file} -o {bam_file} 2>&1 | tee -a {log_file}"
        print("\n", sam_bam_cmd,"\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: sam to bam conversion of mapped file\n")
            handle.write(f"{sam_bam_cmd}\n")
        run = try_except_continue_on_fail(sam_bam_cmd)
        if not run:
            return False

        # sort bam file & calculate depth
        print(f"Running: sorting bam file and calculating depth")
        sort_sam_cmd = f"samtools sort -T {sample_name} {bam_file} -o {sorted_bam_file} " \
                       f"2>&1 | tee -a {log_file}"
        print("\n", sort_sam_cmd, "\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: sorting bam file\n{sort_sam_cmd}\n")
        run = try_except_continue_on_fail(sort_sam_cmd)
        if not run:
            return False
        depth_sam_cmd = f"samtools depth -a {sorted_bam_file} > {depth_file} " \
                       f"2>&1 | tee -a {log_file}"
        print("\n", depth_sam_cmd, "\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: calculating depth\n{depth_sam_cmd}\n")
        run = try_except_continue_on_fail(depth_sam_cmd)
        if not run:
            return False
        positional_depth = {}
        positional_depth_list = []
        with open(depth_file, 'r') as handle:
            for line in csv.reader(handle, dialect="excel-tab"):
                positional_depth[str(line[1])] = int(line[2])
                positional_depth_list.append(int(line[2]))
        if len(positional_depth_list) == 0:
            positional_depth_list.append(0)

        mean_depth = mean(positional_depth_list)
        depth_outfile = pathlib.Path(sample_dir, f"{sample_name}_depth.csv")
        with open(depth_outfile, 'a') as fh:
            fh.write(f"{sample_name}_{ref_name},{mean_depth}\n")

        # index bam file
        print(f"\nRunning: indexing bam file")
        index_bam_cmd = f"samtools index {sorted_bam_file} 2>&1 | tee -a {log_file}"
        print("\n", index_bam_cmd,"\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: indexing bam file\n")
            handle.write(f"{index_bam_cmd}\n")
        run = try_except_continue_on_fail(index_bam_cmd)
        if not run:
            return False

        # convert bam file to a mutli fasta alignment
        print(f"\nRunning: making consensuses sequence from bam to MSA with jvarkit\n")

        sam4web = pathlib.Path(script_dir, "jvarkit", "dist", "sam4weblogo.jar")
        msa_from_bam = f"java -jar {sam4web} -r '{reference_slice}' -o {msa_fasta} " \
                       f"{sorted_bam_file} 2>&1 | tee -a {log_file}"
        print(msa_from_bam)

        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: making consensuses sequence from bam to MSA with jvarkit\n")
            handle.write(f"{msa_from_bam}\n")
        run = try_except_continue_on_fail(msa_from_bam)
        if not run:
            return False

        # convert multi fasta alignment to consensus sequence
        fasta_msa_d = fasta_to_dct(msa_fasta)

        if len(fasta_msa_d) == 0:
            print(f"\nNo MSA made from Bam file\nNo reads may have been mapped\n\n")
            with open(log_file, 'a') as handle:
                handle.write(f"\nNo MSA made from Bam file\nNo reads may have been mapped\n\n")
            empty_file = open(msa_cons, 'w')
            empty_file.close()
            depth_outfile = pathlib.Path(plot_folder, sample_name + '_' + ref_name + "_sequencing_depth.png")
            empty_file = open(depth_outfile, 'w')
            empty_file.close()

        else:
            cons, depth_profile = consensus_maker(fasta_msa_d, positional_depth, min_depth, use_gaps)
            with open(msa_cons, 'w') as handle:
                handle.write(f">{sample_name}_msa\n{cons}\n")

            # write consensus to the sample_viruses file
            with open(sample_viruses_file, 'a') as fh:

                fh.write(f">{sample_name}_{ref_name}_msa\n{cons.replace('-', '')}\n")
                fh.close()

            # plot depth for sample
            depth_list = depth_profile["non_gap"]
            depth_outfile = pathlib.Path(plot_folder, sample_name + '_' + ref_name + "_sequencing_depth.png")
            plot_depth(depth_list, sample_name, depth_outfile)

    # delete single ref files
    for file in pathlib.Path(sample_dir).glob("*fasta*"):
        os.remove(file)

    completed_empty_file = pathlib.Path(sample_dir, sample_name + ".completed")
    empty_file = open(completed_empty_file, 'w')
    empty_file.close()

    print(f"Completed processing sample: {sample_name}\n\n")

    print("Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script runs the read mapping, plotting and consensus generation'
                                                 'for a sample',
                                     formatter_class=Formatter)
    parser.add_argument('-in', '--infile', type=str, default=None, required=True,
                        help='The path and name of the sample fastq file')
    parser.add_argument('-lf', '--log_file', type=str, default=None, required=True,
                        help='The name and path for the logfile')
    parser.add_argument('-rs', '--chosen_ref_file', type=str, default=None, required=True,
                        help='The path of the fasta file with all chosen references')
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="The number of threads to use", required=False)
    parser.add_argument("-d", "--min_depth", type=int, default=100, help="The minimum coverage to call a position in "
                                                                         "the MSA to consensus", required=False)
    parser.add_argument("-ug", "--use_gaps", default=False, action="store_true",
                        help="use gap characters when making the consensus sequences", required=False)

    args = parser.parse_args()
    infile = args.infile
    log_file = args.log_file
    chosen_ref_file = args.chosen_ref_file
    threads = args.threads
    min_depth = args.min_depth
    use_gaps = args.use_gaps

    main(infile, log_file, chosen_ref_file,
         threads, min_depth, use_gaps)
