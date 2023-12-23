import argparse
import pathlib
import os
import subprocess
import csv
from statistics import mean
from src.misc_functions import file_len
from src.misc_functions import try_except_continue_on_fail
from src.misc_functions import consensus_maker
from src.misc_functions import fasta_to_dct
from src.misc_functions import plot_depth
from src.misc_functions import create_coverage_mask

__author__ = 'Philippe Selhorst & Colin Anthony'


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
    plot_folder = pathlib.Path(project_dir, "seq_depth_plots")
    raw_sample_name = sample_name.replace('.no_host', '')
    raw_sample_fastq = pathlib.Path(project_dir, 'raw_samples', f'{raw_sample_name}.fastq')
    sample_viruses_file = pathlib.Path(seq_folder, raw_sample_name + "_viruses.fasta")
    if os.path.isfile(sample_viruses_file):
        os.unlink(sample_viruses_file)
    sam_outfile = pathlib.Path(sample_dir, sample_name + ".ref.sam")
    bam_outfile = pathlib.Path(sample_dir, sample_name + ".ref.sorted.bam")

    # run minimap2 on multi_ref file
    print(f"\nRunning: minimap2 read mapping")
    minimap2_cmd = f"minimap2 --secondary=no -a -Y -t {threads} -x map-ont {chosen_ref_file} " \
                   f"{infile} -o {sam_outfile} 2>&1 | tee -a {log_file}"
    print("\n", minimap2_cmd, "\n")
    with open(log_file, "a") as handle:
        handle.write(f"\nRunning: minimap read mapping\n")
        handle.write(f"{minimap2_cmd}\n")
    run = try_except_continue_on_fail(minimap2_cmd)
    if not run:
        return False

    # convert sam to sorted and indexed bam
    print(f"\nCreating sorted & indexed bam file")
    samtools_cmd = f"samtools view {sam_outfile} -bS -F 2048 | samtools sort - -o {bam_outfile} 2>&1 | tee -a {log_file}\nsamtools index {bam_outfile} 2>&1 | tee -a {log_file}"
    print("\n", samtools_cmd, "\n")
    with open(log_file, "a") as handle:
        handle.write(f"\nRunning: creating sorted & indexed bam file\n")
        handle.write(f"{samtools_cmd}\n")
    run = try_except_continue_on_fail(samtools_cmd)
    if not run:
        return False

    # initialize the fasta files with viruses per sample

    reference_d = fasta_to_dct(chosen_ref_file)
    for k,v in reference_d.items():
        ref_name = k[0:-7]
        ref_seq = v.replace('-', '').lower()
        with open(sample_viruses_file,'a') as handle:
            handle.write(f">{ref_name}\n{ref_seq}\n")
            handle.close()
    depth_outfile1 = pathlib.Path(sample_dir, sample_name + ".depth.csv")
    basecount_file = pathlib.Path(sample_dir, sample_name + ".basecount.csv")
    with open(depth_outfile1, 'a') as fh:
        fh.write(f"sample_name,ref_name,mean_depth,total_reads,virus_reads,percentage\n")


    # iterate mapping over the references
    for k, v in reference_d.items():

        # generate virus dir & set input and output file paths
        ref_name = k[0:-7]
        print(ref_name)
        ref_seq = v.replace('-', '')
        reference_slice = f"{ref_name}:1-{len(ref_seq)}"
        virus_dir = pathlib.Path(sample_dir, ref_name)
        virus_dir.mkdir(mode=0o777, parents=True, exist_ok=True)
        os.chdir(virus_dir)
        single_ref_file = pathlib.Path(virus_dir, ref_name + ".fasta")
        with open(single_ref_file, 'w') as fh:
            fh.write(f">{ref_name}\n{ref_seq}\n")
        depth_file = pathlib.Path(virus_dir, sample_name + f".{ref_name}.depth.tsv")
        reads_file = pathlib.Path(virus_dir, sample_name + f".{ref_name}.reads.txt")
        msa_fasta = pathlib.Path(virus_dir, sample_name + f".{ref_name}.msa_from_bam_file.fasta")
        msa_cons = pathlib.Path(virus_dir, sample_name + f".{ref_name}.msa_consensus.fasta")
        ref_aligned_outfile = pathlib.Path(virus_dir, sample_name + f".{ref_name}.sorted.bam")
        hdf_outfile = pathlib.Path(virus_dir, sample_name + f".{ref_name}.hdf")
        vcf_outfile = pathlib.Path(virus_dir, sample_name + f".{ref_name}.vcf")
        gz_outfile = pathlib.Path(virus_dir, sample_name + f".{ref_name}.vcf.gz")
        vcf_passfile = pathlib.Path(virus_dir, sample_name + f".{ref_name}.pass.vcf")
        gz_passfile = pathlib.Path(virus_dir, sample_name + f".{ref_name}.pass.vcf.gz")
        vcf_failfile = pathlib.Path(virus_dir, sample_name + f".{ref_name}.fail.vcf")
        precon_file = pathlib.Path(virus_dir, sample_name + f".{ref_name}.preconsensus.fasta")
        con_file = pathlib.Path(virus_dir, sample_name + f".{ref_name}.consensus.fasta")

        # extract ref specific alignment from multiref bam file
        with open(log_file, "a") as handle:
            handle.write(f"\n--------{ref_name}--------\n")
        print(f"\nRunning: extracting {ref_name} alignment")
        extract_cmd = f"samtools view -b {bam_outfile} {ref_name} -o - | samtools sort -o {ref_aligned_outfile} 2>&1 | tee -a {log_file}\nsamtools index {ref_aligned_outfile} 2>&1 | tee -a {log_file}"
        print("\n", extract_cmd, "\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: extracting {ref_name} alignment\n")
            handle.write(f"{extract_cmd}\n")
        run = try_except_continue_on_fail(extract_cmd)
        if not run:
            return False

        # # run read mapping using bwa
        # make_index_cmd = f"bwa index {single_ref_file}"
        # with open(log_file, "a") as handle:
        #     handle.write(f"\n{make_index_cmd}\n")
        #
        # try_except_exit_on_fail(make_index_cmd)
        #
        # print(f"\nrunning: bwa read mapping\n")
        # bwa_cmd = f"bwa mem -t {threads} -x ont2d {single_ref_file} {sample_fastq} -o {sam_file} " \
        #           f"2>&1 | tee -a {log_file}"
        # print("\n", bwa_cmd, "\n")
        # with open(log_file, "a") as handle:
        #     handle.write(f"\nrunning: bwa read mapping\n")
        #     handle.write(f"{bwa_cmd}\n")
        # run = try_except_continue_on_fail(bwa_cmd)
        # if not run:
        #     return False


        # # convert sam to bam
        # print(f"\nRunning: sam to bam conversion of mapped file")
        # sam_bam_cmd = f"samtools view -bS {sam_file} -F 2048 -o {bam_file} 2>&1 | tee -a {log_file}"
        # print("\n", sam_bam_cmd,"\n")
        # with open(log_file, "a") as handle:
        #     handle.write(f"\nRunning: sam to bam conversion of mapped file\n")
        #     handle.write(f"{sam_bam_cmd}\n")
        # run = try_except_continue_on_fail(sam_bam_cmd)
        # if not run:
        #     return False

        # sort bam file & calculate depth
        print(f"Running: calculating depth")
        # sort_sam_cmd = f"samtools sort -T {sample_name} {bam_file} -o {sorted_bam_file} " \
        #                f"2>&1 | tee -a {log_file}"
        # print("\n", sort_sam_cmd, "\n")
        # with open(log_file, "a") as handle:
        #     handle.write(f"\nRunning: sorting bam file\n{sort_sam_cmd}\n")
        # run = try_except_continue_on_fail(sort_sam_cmd)
        # if not run:
        #     return False
        depth_sam_cmd = f"samtools depth -a {ref_aligned_outfile} > {depth_file} " \
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


        #get total number of reads and calculate % virus
        total_reads = file_len(raw_sample_fastq)/4
        sam_view_cmd = f"samtools view -F 0x904 -c {ref_aligned_outfile} -o {reads_file} 2>&1 | tee -a {log_file}"
        print("\n", sam_view_cmd, "\n")
        run = try_except_continue_on_fail(sam_view_cmd)
        if not run:
            return False
        with open(reads_file,'r') as fh:
            virus_reads = int(fh.readline())
        percentage = virus_reads/total_reads*100

        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: calculating % virus reads\n{sam_view_cmd}\n")
            handle.write(f"\nTotal reads = {total_reads}\n Virus reads = {virus_reads} \n % virus reads = {percentage}\n")

        with open(depth_outfile1, 'a') as fh:
            fh.write(f"{sample_name},{ref_name},{mean_depth},{total_reads},{virus_reads},{percentage}\n")

        # get total number of bases and calculate % virus
        awk_cmd = f'awk "NR % 4 == 0" ORS="" {raw_sample_fastq}|wc -m'
        print("\n", awk_cmd, "\n")
        total_basecount = int(subprocess.check_output(awk_cmd, shell=True))

        sam_stats_cmd = f'samtools stats -in {ref_aligned_outfile} | grep "bases mapped (cigar):"| cut -f 3'
        print("\n", sam_stats_cmd, "\n")
        basecount = int(subprocess.check_output(sam_stats_cmd, shell=True))
        base_percentage = basecount/total_basecount*100
        with open(basecount_file,'a') as fh:
            fh.write(f"{sample_name},{ref_name},{total_basecount},{basecount},{base_percentage}\n")

        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: calculating % virus bases\n{sam_view_cmd}\n")
            handle.write(f"\nTotal bases = {total_basecount}\n Virus bases = {basecount} \n % virus bases = {base_percentage}\n")

        # # index bam file
        # print(f"\nRunning: indexing bam file")
        # index_bam_cmd = f"samtools index {ref_aligned_outfile} 2>&1 | tee -a {log_file}"
        # print("\n", index_bam_cmd,"\n")
        # with open(log_file, "a") as handle:
        #     handle.write(f"\nRunning: indexing bam file\n")
        #     handle.write(f"{index_bam_cmd}\n")
        # run = try_except_continue_on_fail(index_bam_cmd)
        # if not run:
        #     return False
        
        # convert bam file to a fasta alignment
        print(f"\nRunning: making consensuses sequence from bam to MSA with jvarkit\n")

        sam4web = pathlib.Path(script_dir, "jvarkit", "dist", "sam4weblogo.jar")
        msa_from_bam = f"java -jar {sam4web} -r '{reference_slice}' -o {msa_fasta} " \
                       f"{ref_aligned_outfile} 2>&1 | tee -a {log_file}"
        print(msa_from_bam)

        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: making consensuses sequence from bam to MSA with jvarkit\n")
            handle.write(f"{msa_from_bam}\n")
        run = try_except_continue_on_fail(msa_from_bam)
        if not run:
            return False

        # convert fasta alignment to msa consensus sequence
        fasta_msa_d = fasta_to_dct(msa_fasta)
        sample_name_short = sample_name.split('.')[0]
        ref_name_short = ref_name.split('_')[0]
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
                handle.write(f">{sample_name_short}_{ref_name_short}_msa\n{cons}\n")

            # write consensus to the sample_viruses file
            # with open(sample_viruses_file, 'a') as fh:
            #     fh.write(f">{sample_name}_{ref_name}_msa\n{cons.replace('-', '')}\n")
            #     fh.close()

            # plot depth for sample
            depth_list = depth_profile["non_gap"]
            y = pow(10, int(len(str(total_basecount))))
            print(y)
            norm_depth_list = [x/total_basecount*y for x in depth_list]
            print(depth_list)
            print(norm_depth_list)
            depth_outfile = pathlib.Path(plot_folder, sample_name + '_' + ref_name + "_sequencing_depth.png")
            plot_depth(depth_list, sample_name_short, depth_outfile, ref_name_short)

        # generate artic consensus sequence
        create_coverage_mask(depth_file, 20)
        artic_cmd=[]
        artic_cmd.append(f"medaka consensus --model r1041_e82_400bps_hac_g615 --threads 2 --chunk_len 800 --chunk_ovlp 400 {ref_aligned_outfile} {hdf_outfile} 2>&1 | tee -a {log_file}")
        artic_cmd.append(f"medaka variant {single_ref_file} {hdf_outfile} {vcf_outfile} 2>&1 | tee -a {log_file}")
        artic_cmd.append(f"bgzip -f {vcf_outfile} 2>&1 | tee -a {log_file}")
        artic_cmd.append(f"tabix -f -p vcf {gz_outfile} 2>&1 | tee -a {log_file}")
        artic_cmd.append(f"longshot -P 0 -F -A --no_haps --bam {ref_aligned_outfile} --ref {single_ref_file} --out {vcf_outfile} --potential_variants {gz_outfile} 2>&1 | tee -a {log_file}")
        artic_cmd.append(f"python ~/metatropics/vcf_filter.py --medaka {vcf_outfile} {vcf_passfile} {vcf_failfile} 2>&1 | tee -a {log_file}")
        artic_cmd.append(f"bgzip -f {vcf_passfile}  2>&1 | tee -a {log_file}\ntabix -f -p vcf {gz_passfile} 2>&1 | tee -a {log_file}")
        artic_cmd.append(f"python ~/metatropics/mask.py {single_ref_file} coverage_mask.txt {vcf_failfile} {precon_file} 2>&1 | tee -a {log_file}")
        artic_cmd.append(f"bcftools consensus -f {precon_file} {gz_passfile} -o {con_file} 2>&1 | tee -a {log_file}")

        print("\n", artic_cmd, "\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: artic consensus generation\n")
        for cmd in artic_cmd:
            with open(log_file, "a") as handle:
                handle.write(f"\n{cmd}\n")
            run = try_except_continue_on_fail(cmd)
            if not run:
                return False

        #rename artic consensus header and write to sample viruses file
        with open(con_file, "r") as handle:
            data = handle.readlines()
            data[0] = f">{sample_name_short}_{ref_name_short}_art\n"
        with open(con_file, "w") as handle:
            handle.writelines(data)
        with open(sample_viruses_file, 'a') as fh:
            fh.writelines(data)


    # # delete single ref files
    # for file in pathlib.Path(sample_dir).glob("*fasta*"):
    #     os.remove(file)

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
