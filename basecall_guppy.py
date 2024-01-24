import argparse
import pathlib
from src.misc_functions import try_except_continue_on_fail
from src.json_converter import json_converter
import os, time
import shutil

__author__ = 'Philippe Selhorst'


class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass

def main(inpath, guppy_path, outpath, basecall_mode, real_time, script_folder,barcodes):
    # force absolute file paths
    inpath = pathlib.Path(inpath).absolute()
    outpath = pathlib.Path(outpath).absolute()
    if pathlib.Path.exists(outpath):
        shutil.rmtree(outpath)
    guppy_path = pathlib.Path(guppy_path).absolute()
    guppy_basecaller = pathlib.Path(guppy_path, "guppy_basecaller")
    cuda_device = "CUDA:0"
    gpu_settings = f"-x 'auto' "
    if real_time:
        home = pathlib.Path.home()
        yaml_path = home / "miniconda3/envs/meta_dev/lib/node_modules/artic-rampart/default_protocol/pipelines/demux_map/config.yaml"
        with open(yaml_path, 'r') as file:
            filedata = file.read()
        if barcodes == 'CUST':
            filedata = filedata.replace('barcode_set: "native"', 'barcode_set: "pcr"')
        else:
            filedata = filedata.replace('barcode_set: "pcr"', 'barcode_set: "native"')
        with open(yaml_path, 'w') as file:
            file.write(filedata)

        projectpath = inpath.parent
        json_converter(projectpath, barcodes)
        basecalling_folder = pathlib.Path(projectpath, "basecalling")
        basecalling_folder.mkdir(mode=0o777, parents=True, exist_ok=True)
        temp_folder = pathlib.Path(projectpath, "temp")
        temp_folder.mkdir(mode=0o777, parents=True, exist_ok=True)
        leftover_folder = pathlib.Path(projectpath, "leftover")
        resume = ""
        passfolder = pathlib.Path(projectpath, "fastq/pass")
        passfolder.mkdir(mode=0o777, parents=True, exist_ok=True)
        rampart_protocol_dir = pathlib.Path(script_folder, f"rampart/protocols")
        rampart_protocol_subdirs = [subdir for subdir in os.listdir(rampart_protocol_dir) if os.path.isdir(os.path.join(rampart_protocol_dir,subdir))]
        subdirname=rampart_protocol_subdirs[0]
        protocolpath= pathlib.Path(rampart_protocol_dir,subdirname)
        rampart_cmd = f"rampart --protocol {protocolpath}"
        print(rampart_cmd)
        try_except_continue_on_fail(f"gnome-terminal -- /bin/sh -c 'export NODE_OPTIONS=--max-old-space-size=16384; {rampart_cmd}; exec bash'")
        try_except_continue_on_fail(f"gnome-terminal -- google-chrome http://localhost:3000/")
        counter = 0
        w = 0
        while w == 0:
            pod5files = sorted(os.listdir(inpath), key=lambda y: os.path.getmtime(os.path.join(inpath, y)))
            firstlength = len(pod5files)
            if firstlength > 10:
                x = 10
                counter += x
            else:
                time.sleep(900)
                pod5files = sorted(os.listdir(inpath), key=lambda y: os.path.getmtime(os.path.join(inpath, y)))
                secondlength = len(pod5files)
                if firstlength < secondlength:
                    x = len(pod5files)
                else:
                    x = len(pod5files)
                    w = 1
                counter += x

            for filename in pod5files[0:x]:
                if not filename.startswith('.'):
                    file = os.path.join(inpath, filename)
                    shutil.move(file, basecalling_folder)

            guppy_basecall_cmd = f"{str(guppy_basecaller)} -i {basecalling_folder} -r -s {outpath} -c {basecall_mode} " \
                                 f"--records_per_fastq 4000 --min_qscore 9 {resume}" \
                                 f"{gpu_settings}"

            run = try_except_continue_on_fail(guppy_basecall_cmd)
            if run:
                print(f"Basecalled {counter} pod5 files")
            else:
                print("Basecalling failed")

            for filename in os.listdir(basecalling_folder):
                file = os.path.join(basecalling_folder, filename)
                shutil.move(file, temp_folder)
            resume = '--resume '

        os.rename(inpath, leftover_folder)
        os.rename(temp_folder, inpath)
        os.rmdir(basecalling_folder)

        return True

    else:
        guppy_basecall_cmd = f"{str(guppy_basecaller)} -i {inpath} -r -s {outpath} -c {basecall_mode} " \
                             f"--compress_fastq --records_per_fastq 4000 --min_qscore 9 " \
                             f"{gpu_settings}"
        print(guppy_basecall_cmd)
        run = try_except_continue_on_fail(guppy_basecall_cmd)

        if run:
            print("Basecalling completed\n")
        else:
            print("Basecalling failed")

        return run

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=Formatter)

    parser.add_argument('-in', '--inpath', type=str, default=None, required=True,
                        help='The path to the pod5 folder')
    parser.add_argument('-p', '--guppy_path', type=str, default=None, required=True,
                        help='The path to guppy exexutable')
    parser.add_argument('-sf', '--script_folder', type=str, default=None, required=True,
                        help='The path to script_folder')
    parser.add_argument('-o', '--outpath', type=str, default=None, required=True,
                        help='The path for the outfile')
    parser.add_argument("-b", "--basecall_mode", default="dna_r10.4.1_e8.2_400bps_5khz_hac.cfg",
                        choices=["dna_r10.4.1_e8.2_400bps_5khz_hac.cfg", "dna_r9.4.1_450bps_hac.cfg"], type=str,
                        help="Specify the basecall model given to guppy", required=False)
    parser.add_argument("-rt", "--real_time", default=False, action="store_true",
                        help="start basecalling pod5 files in batches during sequencing", required=False)
    parser.add_argument("-bc", "--barcodes", type=str, choices=["CUST", "SQK-NBD114-24"], required=True,
                        help="Specify barcodes used for porechop demultiplexing")

    args = parser.parse_args()
    inpath = args.inpath
    guppy_path = args.guppy_path
    outpath = args.outpath
    basecall_mode = args.basecall_mode
    real_time = args.real_time
    script_folder = args.scriptfolder
    barcodes = args.barcodes

    main(inpath, guppy_path, outpath, basecall_mode, real_time, script_folder,barcodes)
