import os, time
import sys
import argparse
from pathlib import Path
import subprocess
import datetime
import pandas as pd
import shutil
import glob
import csv

from basecall_dorado import main as dorado_basecall
from demultiplex_dorado import main as dorado_demultiplex
from src.misc_functions import try_except_continue_on_fail
from src.misc_functions import try_except_exit_on_fail
from src.misc_functions import cat_sample_names_filtered
from src.misc_functions import fasta_to_dct
from src.misc_functions import file_len
from src.misc_functions import remove_gaps_in_fasta

__author__ = 'Philippe Selhorst'

class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass

def main(project_dir, min_len, max_len, low_complex, min_depth, run_step,
         rerun_step_only, basecall_mode, cpu_threads,use_gaps, real_time, host, barcodes, one_end):

    # set the dir paths

    script_dir = Path(__file__).absolute().parent
    project_dir = Path(project_dir).absolute()
    reference_seqs_file = Path(script_dir, "references.fasta")
    remove_gaps_in_fasta(reference_seqs_file)
    print(f"\nProject dir is {project_dir}")
    run_name = project_dir.parts[-1]
    pod5_dir = Path(project_dir, "pod5")
    fastq_dir = Path(project_dir, "fastq")
    # pass_dir = Path(fastq_dir, "pass")
    demultiplexed_dir = Path(project_dir, "demultiplexed")
    dorado_dir = ""
    for folder in Path(script_dir).glob("dorado-*-linux-x64/bin"):
        dorado_dir = folder
    nanoplot_dir = Path(project_dir, "nanoplot")
    all_sample_dir = Path(project_dir, "samples")
    no_host_dir = Path(project_dir, "no_host_samples")
    raw_sample_dir = Path(project_dir, "raw_samples")
    sample_names_file = Path(project_dir, "sample_names.csv")
    # seq_folder = Path(project_dir, "seq_files")
    for file in Path(project_dir).glob("*percentages.csv"):
        os.remove(file)
    percentages_file = Path(project_dir, "virus_percentages.csv")
    demulti_host_file = Path(project_dir, "demulti_host.csv")

    seq_summary_file = ""
    for file in project_dir.glob('sequencing_summary*.txt'):
        if not file.exists():
            sys.exit("Could not find sequencing_summary*.txt in project dir")
        else:
            seq_summary_file = file
    plot_folder = Path(project_dir, "seq_depth_plots")
    if os.path.exists(plot_folder):
        shutil.rmtree(plot_folder)

    time_stamp = str('{:%Y-%m-%d_%Hh%M}'.format(datetime.datetime.now()))
    log_file = Path(project_dir, f"{time_stamp}_{run_name}_log_file.txt")
    log_file_final = Path(project_dir, f"{time_stamp}_{run_name}_log_file_final.txt")

    # log start time and change dir to project dir
    print(f"\n# Starting metagenomics pipeline for project: {run_name}\n")
    with open(log_file, "w") as handle:
        handle.write(f"# Starting metagenomics pipeline for project: {run_name}\n")
    now = datetime.datetime.now()
    date_time = now.strftime("%d/%m/%Y, %H:%M:%S")
    print(f"\nstart time = {date_time}\n")
    with open(log_file, "a") as handle:
        handle.write(f"\nstart time = {date_time}\n")
    os.chdir(project_dir)

    # basecalling
    if run_step == 0:
        answer = input("Are you happy with the reference file in the metatropics folder (y/n)?")
        if answer == 'n':
            sys.exit()

        if not sample_names_file.exists():
            sys.exit("Could not find sample_names.csv in project dir")
        print(f"\n________________\n\nRunning: basecalling\n________________\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: basecalling\n")
        run = dorado_basecall(pod5_dir, dorado_dir, fastq_dir, basecall_mode, real_time, script_dir,barcodes)
        # faildir = Path(fastq_dir, "fail")
        # shutil.rmtree(faildir)
        if run and not rerun_step_only:
            run_step = 1
        elif run and rerun_step_only:
            sys.exit("Run step basecalling completed, exiting")
        else:
            sys.exit("Basecalling failed")

    if run_step == 1:
        if not sample_names_file.exists():
            sys.exit("Could not find sample_names.csv in project dir")
        print(f"\n________________\n\nRunning: demultiplexing________________\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: demultiplexing\n")
        if nanoplot_dir.exists():
            os.remove(nanoplot_dir)
        if not list(fastq_dir.glob("*.fastq*")):
            sys.exit(f"No calls.fastq files found in {str(fastq_dir)}")
        else:
            run = dorado_demultiplex(fastq_dir, dorado_dir, demultiplexed_dir, barcodes, one_end)

        if run and not rerun_step_only:
            run_step = 2
        elif run and rerun_step_only:
            sys.exit("Run step demultiplexing completed, exiting")
        else:
            sys.exit("Demultiplexing failed")

    # length filtering and primer trimming allowing for multiple fastqs from multiple exp per barcode
    if run_step == 2:
        if not sample_names_file.exists():
            sys.exit("Could not find sample_names.csv in project dir")

        print("\n________________\n\nRunning: nanoplot, filtering, trimming, and renaming\n________________\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: nanoplot, filtering, trimming, and renaming\n")

        # do Nanoplot
        nanoplot = ""

        if nanoplot_dir.exists():
            answer = input("Previous nanoplot files exists, overwrite (y/n)?")
            if answer == 'y':
                nanoplot = 'yes'
            else:
                nanoplot = 'no'
                print(f"\nKeeping previous NanoPlot output\n")
                with open(log_file, "a") as handle:
                    handle.write(f"\nKeeping previous NanoPlot output\n")
        else:
            nanoplot = 'yes'
        if nanoplot == 'yes':
            print(f"\nPerforming NanoPlot on demultiplexed fastqs\n")
            with open(log_file, "a") as handle:
                handle.write(f"\nPerforming NanoPlot on demultiplexed fastqs\n")
            for folder in project_dir.glob("nanoplot"):
                shutil.rmtree(folder)
            all_nanoplot_files = demultiplexed_dir.glob("*_barcode*")
            for barcode_fastq in all_nanoplot_files:
                sample_name = Path(barcode_fastq).stem
                nanoplotoutputdir = Path(project_dir, f"nanoplot/{sample_name}")
                nanoplotcmd = f"NanoPlot --fastq {barcode_fastq} -o {nanoplotoutputdir}"
                print(nanoplotcmd)
                subprocess.call(nanoplotcmd, shell=True)

        raw_sample_dir.mkdir(mode=0o777, parents=True, exist_ok=True)

        classified_reads = 0

        unclassified_file  = list(Path(demultiplexed_dir).glob("*unclassified.fastq"))[0]
        unclassified_reads = file_len(unclassified_file) / 4
        pre_existing_files = list(demultiplexed_dir.glob("*_bad.*"))
        if pre_existing_files:
            answer = input("Previous filtered files exists, overwrite (y/n)?")
            if answer == 'n':
                sys.exit("\nKeeping filtered files, exiting\n")
            else:
                for file in pre_existing_files:
                    os.remove(file)

        print(f"\nFiltering and trimming with following parameters\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nFiltering and trimming with following parameters\n")
        left_trim = ''
        right_trim = ''
        lowcomp = ''

        if barcodes == "CUST" or barcodes == "SQK-RPB114-24":
            print(f"min_length = {min_len},max_length = {max_len}, trimmed 0bp")
            with open(log_file, "a") as handle:
                handle.write(f"min_length = {min_len},max_length = {max_len}, trimmed 0bp")
        else:
            print(f"min_length = {min_len},max_length = {max_len}, trimmed 18bp solA")
            with open(log_file, "a") as handle:
                handle.write(f"min_length = {min_len},max_length = {max_len}, trimmed 18bp solA")
            left_trim = '-trim_left 18'
            right_trim = '-trim_right 18'
        if low_complex == '-lc':
            lowcomp = '-lc_dust 0.1 -derep'
            print (f"0.1 low complexity and replicates filter applied")
            with open(log_file, "a") as handle:
                handle.write(f"0.1 low complexity and replicates filter applied")

        for file in demultiplexed_dir.glob("*_barcode*"):

            barcode_number = file.parts[-1].split('.')[0].split('_')[-1]
            classified_reads += (file_len(file) / 4)
            out_good = Path(demultiplexed_dir, f"{barcode_number}_good.fastq")
            out_bad = Path(demultiplexed_dir, f"{barcode_number}_bad.fastq")
            prinseq_cmd = f"source $(conda info --base)/etc/profile.d/conda.sh && conda activate prinseq-plus-plus && prinseq++ -fastq {file} -threads {cpu_threads} {lowcomp} {left_trim} {right_trim} -min_len {min_len} -max_len {max_len} -VERBOSE=1 -out_bad {out_bad} -out_good {out_good} 2>&1 | tee -a {log_file} && conda activate meta"
            print(f"\nFiltering and trimming {file}\n {prinseq_cmd}\n")
            subprocess.call(prinseq_cmd, shell=True, executable="/bin/bash")
            if not out_good:
                print(f"No sequences in file after filtering and primer trimming for {file}\n")

        percentage_unclassified = unclassified_reads/(classified_reads+unclassified_reads)*100
        print(f"\nPercentage unclassified reads is {percentage_unclassified}\n")
        with open(demulti_host_file, 'w') as fh:
            fh.write(f"percentage_unclassified,{percentage_unclassified}\n")

        # do rename
        print(f"\nRenaming and concatenating\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nRenaming and concatenating\n")
        sample_names_df = pd.read_csv(sample_names_file, sep=None, keep_default_na=False, na_values=['NA'],
                                      engine="python")
        sample_names_df['barcode_1'] = sample_names_df['barcode_1'].apply(lambda x: cat_sample_names_filtered(x))
        sample_names_df['barcode_2'] = sample_names_df['barcode_2'].apply(lambda x: cat_sample_names_filtered(x))
        sample_names_dict = sample_names_df.set_index('sample_name').T.to_dict(orient='list')

        for sample_name, [barcode_1, barcode_2] in sample_names_dict.items():

            # sample_dir = Path(all_sample_dir, sample_name)
            # if not sample_dir.exists():
            #     Path(sample_dir).mkdir(mode=0o777, parents=True, exist_ok=True)

            barcode_1_file = Path(demultiplexed_dir, barcode_1)
            # allow for case where only one barcode was specified per sample.
            if barcode_2 == " ":
                barcode_2_file = ""
            else:
                barcode_2_file = Path(demultiplexed_dir, barcode_2)
            cat_outfile = Path(raw_sample_dir, f"{sample_name}.fastq")
            cat_cmd = f"cat {str(barcode_1_file)} {str(barcode_2_file)} > {cat_outfile} 2>&1 | tee -a {log_file}"
            print(cat_cmd)
            with open(log_file, "a") as handle:
                handle.write(f"\n{cat_cmd}\n")
            run = try_except_continue_on_fail(cat_cmd)
            if not run:
                print("Missing one or more demultiplexed files for this sample")
                with open(log_file, "a") as handle:
                    handle.write("\nMissing one or more demultiplexed files for this sample\n")
                continue

        filtered_files = list(demultiplexed_dir.glob("*good.fastq"))
        for file in filtered_files:
            file.unlink()

        if not rerun_step_only:
            run_step = 3
        elif rerun_step_only:
            sys.exit("Nanoplot, filtering, trimming, and renaming completed, exiting")
        else:
            sys.exit("Nanoplot, filtering, trimming, and renaming failed")


    if run_step == 3:
        print("\n________________\n\nRunning: host removal using minimap2\n________________\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: host removal using minimap2\n")
        pre_existing_files = list(raw_sample_dir.glob("*.fastq"))
        if not pre_existing_files:
            sys.exit("No files found in raw sample folder\n")

        if host != '':
            no_host_dir.mkdir(mode=0o777, parents=True, exist_ok=True)
            host_dir = Path(script_dir, "host_genomes", host)
            host_name = list(host_dir.glob("*.fasta"))[0]
            with open(demulti_host_file, 'a') as fh:
                fh.write(f"host,{host_name}\n")

            rib_ref = list(host_dir.glob("18S.fa"))[0]
            print(rib_ref)
            print(f'Host genome to remove is {host_name}')

            with open(demulti_host_file, 'a') as fh:
                fh.write(f"sample_name,total_reads,percentage_18S\n")
            for file in pre_existing_files:
                total_reads = file_len(file) / 4
                sample_name = file.stem
                rib_ref_outfile = Path(raw_sample_dir, f"{sample_name}.18S.fastq")
                minimap_cmd = f"minimap2 --secondary=no -a -Y -t 15 -x map-ont {rib_ref} {file} | samtools view -bF 2308 - | samtools fastq - > {rib_ref_outfile}"
                print(minimap_cmd)
                run = try_except_continue_on_fail(minimap_cmd)
                if not run:
                    print("18S check failed")
                    with open(log_file, "a") as handle:
                        handle.write("\n18S check failed\n")
                    continue
                rib_ref_reads = file_len(rib_ref_outfile)/4
                percentage_rib_ref = (int(rib_ref_reads) / int(total_reads)) * 100
                with open(demulti_host_file, 'a') as fh:
                    fh.write(f"{sample_name},{total_reads},{percentage_rib_ref}\n")
                rib_ref_outfile.unlink()
            with open(demulti_host_file, 'a') as fh:
                fh.write(f"sample_name,total_reads,percentage_host\n")
            for file in pre_existing_files:
                total_reads = file_len(file) / 4
                sample_name = file.stem
                unmapped_outfile = Path(no_host_dir, f"{sample_name}.no_host.fastq")
                minimap_cmd = f"minimap2 --secondary=no -a -Y -t 15 -x map-ont {host_name} {file} | samtools view -f4 - | samtools fastq - > {unmapped_outfile}"
                print(minimap_cmd)
                with open(log_file, "a") as handle:
                    handle.write(f"\n{minimap_cmd}\n")
                run = try_except_continue_on_fail(minimap_cmd)
                if not run:
                    print("Host removal failed")
                    with open(log_file, "a") as handle:
                        handle.write("\nHost removal failed\n")
                    continue
                else:
                    unmapped_reads = file_len(unmapped_outfile)/4
                    percentage_host = (1-(int(unmapped_reads)/int(total_reads)))*100
                    with open(demulti_host_file, 'a') as fh:
                        fh.write(f"{sample_name},{total_reads},{percentage_host}\n")
        else:
            print(f'No host genome to remove')
            with open(log_file, "a") as handle:
                handle.write(f"\nNo host genome to remove\n")


        if not rerun_step_only:
            run_step = 4
        elif rerun_step_only:
            sys.exit("Host removal using minimap2 completed, exiting")
        else:
            sys.exit("Host removal using minimap2 failed")

    # Reference-based assembly
    if run_step == 4:

        if not all_sample_dir.exists():
            all_sample_dir.mkdir(mode=0o777, parents=True, exist_ok=True)
        os.chdir(all_sample_dir)

        pre_existing_files = list(raw_sample_dir.glob("*.fastq"))
        if not pre_existing_files:
            sys.exit("No files found in raw sample folder, exiting\n")
        if host != '':
            pre_existing_files = list(no_host_dir.glob("*.fastq"))
            if not pre_existing_files:
                sys.exit("No files found in no host sample folder, exiting\n")


        # delete pre existing files in project dir
        # for file in Path(project_dir).glob("*.fasta"):
        #     os.remove(file)

        for file in Path(project_dir).glob("*.txt"):
            if "msa" in str(file):
                os.remove(file)

        # delete pre existing virus folders
        for folder in glob.glob("*/*/"):
            shutil.rmtree(folder)

        # delete pre existing files in sample folders
        for file in all_sample_dir.glob("*/*.*"):
            os.remove(file)


        for file in raw_sample_dir.glob("*.fastq"):
            sample_name = file.stem
            sample_dir = Path(all_sample_dir, sample_name)
            if not sample_dir.exists():
                Path(sample_dir).mkdir(mode=0o777, parents=True, exist_ok=True)


        print("\n________________\n\nRunning: reference-based assembly\n________________\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: reference-based assembly\n")

        # get number of samples and threads
        number_samples = (len(list(raw_sample_dir.glob('*.fastq'))))
        print("number of samples=" + str(number_samples))
        max_threads = cpu_threads
        used_threads = 0
        msa_threads = 2

        log_file_msa_temp = Path(project_dir, f"{time_stamp}_{run_name}_log_file_msa_temp.txt")
        log_file_msa = Path(project_dir, f"{time_stamp}_{run_name}_log_file_msa.txt")

        print(f"min_depth = {min_depth}")
        with open(log_file_msa_temp, "a") as handle:
            handle.write(f"\nmin_depth = {min_depth}\n")

        if host != '':
            all_sample_files = Path(no_host_dir).glob("*.fastq")
        else:
            all_sample_files = Path(raw_sample_dir).glob("*.fastq")
        sample_no = 0
        for sample_fastq in all_sample_files:

            sample_no += 1

            # get fastq path and name
            sample_dir = Path(sample_fastq).parent
            sample_name = Path(sample_fastq).stem
            log_file_msa_sample = Path(project_dir, f"{time_stamp}_{sample_name}_log_file_msa_sample.txt")
            os.chdir(sample_dir)

            # check free threads
            finished_threads = len(list(Path(all_sample_dir).glob("*/*.completed")))*msa_threads
            free_threads = max_threads + finished_threads - used_threads
            print("\nfree_threads = " + str(free_threads))
            print('\n' + 'waiting for free threads')
            while free_threads < msa_threads:
                time.sleep(5)
                finished_threads = len(list(Path(all_sample_dir).glob("*/*.completed")))*msa_threads
                free_threads = max_threads + finished_threads - used_threads
                print(free_threads)

            # check if fastq is present
            file_present = list(sample_dir.glob("*.fastq"))
            if not file_present:
                print(
                    f"\nCould not find the concatenated sample fastq file in sample folder: {sample_fastq}\nskipping sample")
                with open(log_file_msa_sample, "a") as handle:
                    handle.write(
                        f"\nCould not find the concatenated sample fastq file in sample folder: {sample_fastq}\nskipping sample")
                continue
            print(f"\n------->Running majority consensus pipeline for {sample_no} st/nd sample {sample_name} in new window\n")
            with open(log_file_msa_sample, "a") as handle:
                handle.write(
                    f"\n\n------->Running majority consensus pipeline for {sample_no} st/nd sample {sample_name} in new window\n")

            # start majority consensus pipeline in new window
            majority_cmd = f"python ~/metatropics/msa_consensus.py -in {sample_fastq} -lf {log_file_msa_sample} " \
                           f"-rs {reference_seqs_file} " \
                           f"-t {msa_threads} -d {min_depth} {use_gaps} -b {basecall_mode}"
            print(majority_cmd)
            try_except_continue_on_fail(f"gnome-terminal -- /bin/sh -c 'conda run -n meta {majority_cmd}'")
            used_threads += msa_threads

        # concat all log files
        finished_threads = len(list(Path(all_sample_dir).glob("*/*.completed")))
        while finished_threads < number_samples:
            time.sleep(5)
            finished_threads = len(list(Path(all_sample_dir).glob("*/*.completed")))

        else:
            os.chdir(project_dir)
            loglist_msa = []
            for path in Path(project_dir).glob("*_log_file_msa_sample.txt"):
                loglist_msa.append(str(path))
            sep = " "
            string_msa = sep.join(loglist_msa)
            cat_cmd = f"cat {str(log_file_msa_temp)} {string_msa} > {log_file_msa}"
            try_except_continue_on_fail(cat_cmd)
            cat_cmd = f"cat {str(log_file)} {str(log_file_msa)} > {log_file_final}"
            try_except_continue_on_fail(cat_cmd)
            for path in list(Path(project_dir).glob("*_log_file_msa_sample.txt")):
                os.remove(path)
            os.remove(log_file_msa_temp)
            os.remove(log_file_msa)
            os.remove(log_file)


        #collect & concat all csv files
        os.chdir(all_sample_dir)
        viruslist=['sample_name', '']
        virusdct= fasta_to_dct(reference_seqs_file)
        for virusname, sequence in virusdct.items():
            viruslist.append(virusname[0:-7])

        with open(percentages_file, 'a') as fh:
            fh.write("\n")
            writer = csv.writer(fh)
            writer.writerow(viruslist)

        for csvfile in sorted(Path(all_sample_dir).glob("*/*basecount.csv")):
            opencsv = open(csvfile, 'r')
            csvfile_stem = csvfile.stem.split(".")[0]
            counts = [csvfile_stem, 'base_count']
            percentage = [csvfile_stem, 'base_percent']
            avg_length = [csvfile_stem, 'avg_length']
            for line in csv.reader(opencsv):
                counts.append(line[4])
                percentage.append(line[5])
                avg_length.append(line[2])
            with open(percentages_file, 'a') as fh:
                writer = csv.writer(fh)
                writer.writerow(counts)
                writer.writerow(percentage)
                writer.writerow(avg_length)
            opencsv.close()

        with open(percentages_file, 'a') as fh:
            writer = csv.writer(fh)
            writer.writerow("\n")

        for csvfile in sorted(Path(all_sample_dir).glob("*/*depth.csv")):
            opencsv = open(csvfile, 'r')
            csvfile_stem = csvfile.stem.split(".")[0]
            counts = [csvfile_stem, 'read_count']
            percentage = [csvfile_stem, 'read_percent']
            for line in csv.reader(opencsv):
                if line[0] !="sample_name":
                    counts.append(line[4])
                    percentage.append(line[5])
            with open(percentages_file, 'a') as fh:
                writer = csv.writer(fh)
                writer.writerow(counts)
                writer.writerow(percentage)
            opencsv.close()

        # print("Aligning consensus sequences\n")

        # for seqfile in Path(seq_folder).glob("*.fasta"):
        #     seqfile_name = Path(seqfile).stem
        #     tmp_file = Path(seq_folder, seqfile_name + "_temp_aligned.fasta")
        #     mafft_cmd = f"mafft --thread -1 --auto {str(seqfile)} > {str(tmp_file)}"
        #     print(mafft_cmd)
        #     run = try_except_continue_on_fail(mafft_cmd)
        #     if not run:
        #         print(f"could not align {seqfile}")
        #         sys.exit("exiting")
        #     else:
        #         seqfile.unlink()
        #         os.rename(str(tmp_file), str(seqfile))

    # print end time
    now = datetime.datetime.now()
    date_time = now.strftime("%d/%m/%Y, %H:%M:%S")
    print(f"\nend time = {date_time}\n\n")
    with open(log_file_final, "a") as handle:
        handle.write(f"\nend time = {date_time}\n\n")

    print("Sample processing completed\n")
    with open(log_file_final, "a") as handle:
        handle.write(f"\nSample processing completed\n\n")
    for file in Path(all_sample_dir).glob('*/*.completed'):
        os.remove(file)

    # compress pod5 files
    os.chdir(project_dir)
    targzpath = Path(project_dir.parent, run_name + ".tar")
    pod5_dir_name = pod5_dir.parts[-1]
    seq_summary_file_name = Path(seq_summary_file).name
    tarcmd = f"tar -cf {targzpath} {pod5_dir_name} {seq_summary_file_name}"
    print(tarcmd)
    try_except_exit_on_fail(tarcmd)
    zipcmd = f"pigz -7 -p 16 {targzpath}"
    try_except_exit_on_fail(zipcmd)


    with open(log_file_final, "a") as handle:
        handle.write(f"\n{tarcmd}\n\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process raw nanopore reads to fasta consensus sequences v2.0",
                                     formatter_class=Formatter)
    parser.add_argument("-in", "--project_dir", default=argparse.SUPPRESS, type=str,
                        help="The path to the directory containing the 'pod5' and 'fastq' dirs ", required=True)
    parser.add_argument("-mi", "--min_len", type=int, default=150,
                        help="The minimum read length allowed", required=False)
    parser.add_argument("-ma", "--max_len", type=int, default=1000000,
                        help="The maximum read length allowed", required=False)
    parser.add_argument("-lc", "--low_complex", default='', action="store_const", const='-lc',
                        help="use prinseq low complexity filter", required=False)
    parser.add_argument("-d", "--min_depth", type=int, default=100, help="The minimum coverage to call a position in the MSA to consensus", required=False)
    parser.add_argument("--run_step", default=0, type=int, required=False,
                        help="Run the pipeline starting at this step:\n"
                             "--run_step 0 = basecall reads with dorado\n"
                             "--run_step 1 = demultiplex reads with dorado\n"
                             "--run_step 2 = concatenate, filtering, trimming, rename, combine barcodes, nanoplot\n"
                             "--run_step 3 = remove host reads from sample files\n"
                             "--run_step 4 = run reference-based viral genome assembly on each sample\n")
    parser.add_argument("--run_step_only", default=False, action="store_true",
                        help="Only run the step specified in 'run_step'", required=False)
    parser.add_argument("-b", "--basecall_mode", default="dna_r10.4.1_e8.2_400bps_5khz_hac.cfg", choices=["dna_r10.4.1_e8.2_400bps_5khz_hac.cfg", "dna_r9.4.1_450bps_hac.cfg"], type=str,
                        help="Specify the basecall model given to dorado", required=False)
    parser.add_argument("-c", "--cpu_threads", type=int, default=16, choices=range(0, 21),
                        help="The number of cpu threads to use", required=False)
    parser.add_argument("-ug", "--use_gaps", default='', action="store_const", const='-ug',
                        help="use gap characters when making the consensus sequences", required=False)
    parser.add_argument("-rt", "--real_time", default=False, action="store_true",
                        help="start basecalling pod5 files in batches during sequencing", required=False)
    parser.add_argument("-ho", "--host", default='', type=str, choices=["homo_sapiens","culex","mastomys_natalensis", "mus_musculus", "bos_taurus"], required=False,
                        help="name of host species to remove")
    parser.add_argument("-bc", "--barcodes", type=str, choices=["CUST","SQK-NBD114-24", "SQK-RPB114-24"], required=True,
                        help="Specify barcodes used for demultiplexing, if NBC, 27bp are trimmed from both ends of each read after demultiplexing")
    parser.add_argument("-oe", "--one_end", default=False, action="store_true", required=False,
                        help="use reads if they have barcode on only one end, this increases the amount of data yet increases probability of misclassification")

    args = parser.parse_args()

    project_dir = args.project_dir
    min_len = args.min_len
    max_len = args.max_len
    low_complex= args.low_complex
    min_depth = args.min_depth
    run_step = args.run_step
    run_step_only = args.run_step_only
    basecall_mode = args.basecall_mode
    cpu_threads = args.cpu_threads
    use_gaps = args.use_gaps
    real_time = args.real_time
    host = args.host
    barcodes = args.barcodes
    one_end = args.one_end

    main(project_dir, min_len, max_len, low_complex, min_depth, run_step,
         run_step_only, basecall_mode, cpu_threads, use_gaps, real_time, host, barcodes, one_end)

