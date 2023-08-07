import os, time
import sys
import argparse
from pathlib import Path
import subprocess
import datetime
import pandas as pd
import shutil
import glob
from basecall_guppy import main as gupppy_basecall
from demultiplex_guppy import main as guppy_demultiplex
from src.misc_functions import try_except_continue_on_fail
from src.misc_functions import try_except_exit_on_fail
from src.misc_functions import cat_sample_names
from src.misc_functions import filter_length_trim_seq

__author__ = 'Philippe Selhorst'

class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass

def main(project_dir, reference, ref_start, ref_end, min_len, max_len, min_depth, run_step,
         rerun_step_only, basecall_mode, cpu_threads,use_gaps,
         guppy_dir, real_time):

    # set the dir paths
    script_dir = Path(__file__).absolute().parent
    project_dir = Path(project_dir).absolute()
    reference_seqs_file = Path(script_dir, "references.fasta")
    print(f"\nProject dir is {project_dir}")
    run_name = project_dir.parts[-1]
    fast5_dir = Path(project_dir, "fast5")
    fastq_dir = Path(project_dir, "fastq")
    demultiplexed_dir = Path(project_dir, "demultiplexed")
    all_sample_dir = Path(project_dir, "samples")
    sample_names_file = Path(project_dir, "sample_names.csv")
    if not sample_names_file:
        sys.exit("Could not find sample_names.csv in project dir")
    seq_summary_file = ""
    for file in project_dir.glob('sequencing_summary*.txt'):
        if not file:
            sys.exit("Could not find sequencing_summary*.txt in project dir")
        else:
            seq_summary_file = file
    plot_folder = Path(project_dir, "seq_depth_plots")
    if os.path.exists(plot_folder):
        shutil.rmtree(plot_folder)
    plot_folder.mkdir(mode=0o777, parents=True, exist_ok=True)

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
        print(f"\n________________\n\nRunning: basecalling\n________________\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: basecalling\n")
        run = gupppy_basecall(fast5_dir, guppy_dir, fastq_dir, basecall_mode, real_time, script_dir)
        faildir = Path(fastq_dir, "fail")
        shutil.rmtree(faildir)
        if run and not rerun_step_only:
            run_step = 1
        elif run and rerun_step_only:
            sys.exit("Run step basecalling completed, exiting")
        else:
            sys.exit("Basecalling failed")

    if run_step == 1:
        print(f"\n________________\n\nRunning: demultiplexing________________\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: demultiplexing")
        if not list(fastq_dir.glob("*.fastq*")):
            fastq_dir = Path(fastq_dir, "pass")
            if not list(fastq_dir.glob("*.fastq*")):
                sys.exit(f"No fastq files found in {str(fastq_dir)} or {str(fastq_dir.parent)}")
        run = guppy_demultiplex(fastq_dir, guppy_path, demultiplexed_dir)
        if run and not rerun_step_only:
            run_step = 2
        elif run and rerun_step_only:
            sys.exit("Run step demultiplexing completed, exiting")
        else:
            sys.exit("Demultiplexing failed")

    # length filtering and primer trimming allowing for multiple fastqs from multiple exp per barcode
    if run_step == 2:
        print("\n________________\n\nRunning: length filtering & primer trim\n________________\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nRunning: length filtering & primer trim\n")
        pre_existing_files = list(demultiplexed_dir.glob("*.fastq"))
        if pre_existing_files:
            print("Found existing files in demultiplex folder.\nThese files will be deleted\n")
            for file in pre_existing_files:
                file.unlink()
        for folder in demultiplexed_dir.glob("barcode*"):
            search = list(Path(folder).glob("*.fastq"))
            if not search:
                print(f"No files in folder\nSkipping folder: {folder}\n")
                continue
            if len(search) > 1:
                print(f"Length filtering and trimming {folder}")
                barcode_number = Path(search[0]).parent.parts[-1]
                concat_outfile = f"cat_barcode_{barcode_number}.fastq"
                cat_cmd = f"cat "
                for file in search:
                    cat_cmd += f"{str(file)} "
                cat_cmd += f" > {concat_outfile}"
                try_except_exit_on_fail(cat_cmd)
                new_name = Path(demultiplexed_dir, f"{run_name}_{barcode_number}.fastq")
                filtered_file = filter_length_trim_seq(concat_outfile, new_name, max_len, min_len, 27, 27)
                os.unlink(str(concat_outfile))
                if not filtered_file:
                    print(f"No sequences in file after length filtering and primer trimming for {concat_outfile}\n")
            else:
                print(f"Length filtering and primer trimming {folder}")
                file = Path(search[0])
                barcode_number = file.parent.parts[-1]
                new_name = Path(demultiplexed_dir, f"{run_name}_{barcode_number}.fastq")
                filtered_file = filter_length_trim_seq(file, new_name, max_len, min_len, 27, 27)
                if not filtered_file:
                    print(f"No sequences in file after length filtering and primer trimming for {file}\n")
        if not rerun_step_only:
            run_step = 3
        elif rerun_step_only:
            sys.exit("\nFiltering, trimming, and renaming of demultiplexed files completed, exiting\n")
        else:
            sys.exit("\nFiltering, trimming and renaming of demultiplexed files failed\n")

    # concatenate multiple barcodes per sample and rename according to sample file
    if run_step == 3:
        print("\n________________\n\nRunning: concatenate & rename\n________________\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nCollecting demultiplexed files into sample.fastq files based on specified sample "
                         f"barcode combinations\n")

        sample_names_df = pd.read_csv(sample_names_file, sep=None, keep_default_na=False, na_values=['NA'], engine="python")
        sample_names_df['barcode_1'] = sample_names_df['barcode_1'].apply(lambda x: cat_sample_names(x, run_name))
        sample_names_df['barcode_2'] = sample_names_df['barcode_2'].apply(lambda x: cat_sample_names(x, run_name))
        sample_names_dict = sample_names_df.set_index('sample_name').T.to_dict(orient='list')

        for sample_name, [barcode_1, barcode_2] in sample_names_dict.items():
            sample_dir = Path(all_sample_dir, sample_name)
            if not sample_dir.exists():
                Path(sample_dir).mkdir(mode=0o777, parents=True, exist_ok=True)
            barcode_1_file = Path(demultiplexed_dir, barcode_1)
            # allow for case where only one barcode was specified per sample.
            if barcode_2 == " ":
                barcode_2_file = ""
            else:
                barcode_2_file = Path(demultiplexed_dir, barcode_2)
            cat_outfile = Path(sample_dir, f"{sample_name}.fastq")
            cat_cmd = f"cat {str(barcode_1_file)} {str(barcode_2_file)} > {cat_outfile}"
            print(cat_cmd)
            with open(log_file, "a") as handle:
                handle.write(f"\n{cat_cmd}\n")
            run = try_except_continue_on_fail(cat_cmd)
            if not run:
                print("Missing one or more demultiplexed files for this sample")
                with open(log_file, "a") as handle:
                    handle.write("\nMissing one or more demultiplexed files for this sample\n")
                continue
        if not rerun_step_only:
            run_step = 4
        else:
            sys.exit("Run step only completed, exiting")

    # Reference-based assembly
    if run_step == 4:

        os.chdir(all_sample_dir)

        # delete pre existing files in project dir
        for file in Path(project_dir).glob("*.fasta"):
            os.remove(file)
        for file in Path(project_dir).glob("*.txt"):
            if not "sequencing_summary" in str(file):
                os.remove(file)

        print("\n________________\n\nRunning: reference-based assembly\n________________\n")
        with open(log_file, "a") as handle:
            handle.write(f"\nStarting reference-based assembly\n")

        # get number of samples and threads
        number_samples = (len(list(all_sample_dir.glob('*/*.fastq'))))
        print("number of samples=" + str(number_samples))
        max_threads = cpu_threads
        used_threads = 0
        msa_threads = 2

        # delete pre existing virus folders
        for folder in glob.glob("*/*/"):
            shutil.rmtree(folder)
        # delete pre existing files
        for file in Path(all_sample_dir).glob("*/*.*"):
            if not str(file).endswith(".fastq"):
                os.remove(file)

        log_file_msa_temp = Path(project_dir, f"{time_stamp}_{run_name}_log_file_msa_temp.txt")
        log_file_msa = Path(project_dir, f"{time_stamp}_{run_name}_log_file_msa.txt")

        print(f"min_depth = {min_depth}")
        with open(log_file_msa_temp, "a") as handle:
            handle.write(f"\nmin_depth = {min_depth}\n")

        all_sample_files = Path(all_sample_dir).glob("*/*.fastq")
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

            # do Nanoplot
            nanoplotoutputdir= Path(project_dir, f"nanoplot/{sample_name}")
            nanoplotcmd = f"NanoPlot --fastq {sample_fastq} -o {nanoplotoutputdir}"
            print(nanoplotcmd)
            subprocess.call(nanoplotcmd, shell=True)


            # start majority consensus pipeline in new window
            majority_cmd = f"python ~/metatropics/msa_consensus.py -in {sample_fastq} -lf {log_file_msa_sample} " \
                           f"-rs {reference_seqs_file} " \
                           f"-t {msa_threads} -d {min_depth} {use_gaps}"
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

    # compress fast5 files
    os.chdir(project_dir)
    targzpath = Path(project_dir.parent, run_name + ".tar")
    fast5_dir_name = fast5_dir.parts[-1]
    seq_summary_file_name = Path(seq_summary_file).name
    tarcmd = f"tar -cf {targzpath} {fast5_dir_name} {seq_summary_file_name}"
    print(tarcmd)
    try_except_exit_on_fail(tarcmd)
    zipcmd = f"pigz -7 -p 16 {targzpath}"
    try_except_exit_on_fail(zipcmd)



    with open(log_file_final, "a") as handle:
        handle.write(f"\n{tarcmd}\n\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process raw nanopore reads to fasta consensus sequences",
                                     formatter_class=Formatter)

    parser.add_argument("-in", "--project_dir", default=argparse.SUPPRESS, type=str,
                        help="The path to the directory containing the 'fast5' and 'fastq' dirs ", required=True)
    parser.add_argument("-r", "--reference", type=str, help="The reference genome and primer scheme to use",
                        choices=["ChikAsian_V1_400", "ChikECSA_V1_800", "ZikaAsian_V1_400", "SARS2_V1_800", "SARS2_V1_400", "DENV1_V1_400", "DENV2_V1_400"], required=False)
    parser.add_argument("-rs", "--reference_start", default=1, type=int,
                        help="The start coordinate of the reference sequence for read mapping", required=False)
    parser.add_argument("-re", "--reference_end", default=False, type=int,
                        help="The end coordinate of the reference sequence for read mapping. Default = full length",
                        required=False)
    parser.add_argument("-mi", "--min_len", type=int, default=150,
                        help="The minimum read length allowed:\n = 300 for 400bp amplicon design"
                                                             "\n = 700 for 800bp amplicon design", required=False)
    parser.add_argument("-ma", "--max_len", type=int, default=1000000,
                        help="The maximum read length allowed:\n = 500 for 400bp amplicon design"
                             "                                \n = 900 for 800bp amplicon design", required=False)
    parser.add_argument("-d", "--min_depth", type=int, default=100, help="The minimum coverage to call a position in "
                                                                         "the MSA to consensus", required=False)
    parser.add_argument("--run_step", default=0, type=int, required=False,
                        help="Run the pipeline starting at this step:\n"
                             "--run_step 0 = basecall reads with Guppy\n"
                             "--run_step 1 = demultiplex with Guppy\n"
                             "--run_step 2 = size filer and rename demultiplexed fastq file\n"
                             "--run_step 3 = concatenate demultiplexed files into sample files\n"
                             "--run_step 4 = run reference-based viral genome assembly on each sample\n")
    parser.add_argument("--run_step_only", default=False, action="store_true",
                        help="Only run the step specified in 'run_step'", required=False)
    parser.add_argument("-b", "--basecall_mode", default=1, choices=[0, 1], type=int,
                        help="0 = basecall in fast mode\n"
                             "1 = basecall in high accuracy mode\n", required=False)
    parser.add_argument("-c", "--cpu_threads", type=int, default=14, choices=range(0, 16),
                        help="The number of cpu threads to use", required=False)
    parser.add_argument("-ug", "--use_gaps", default='', action="store_const", const='-ug',
                        help="use gap characters when making the consensus sequences", required=False)
    parser.add_argument("-p", "--guppy_path", default=argparse.SUPPRESS, type=str,
                        help="The path to the guppy executables eg: '.../ont-guppy/bin/'", required=True)
    parser.add_argument("-rt", "--real_time", default=False, action="store_true",
                        help="start basecalling fast5 files in batches during sequencing", required=False)

    args = parser.parse_args()

    project_dir = args.project_dir
    reference = args.reference
    reference_start = args.reference_start
    reference_end = args.reference_end
    min_len = args.min_len
    max_len = args.max_len
    min_depth = args.min_depth
    run_step = args.run_step
    run_step_only = args.run_step_only
    basecall_mode = args.basecall_mode
    cpu_threads = args.cpu_threads
    use_gaps = args.use_gaps
    guppy_path = args.guppy_path
    real_time = args.real_time

    main(project_dir, reference, reference_start, reference_end, min_len, max_len, min_depth, run_step,
         run_step_only, basecall_mode, cpu_threads, use_gaps, guppy_path, real_time)