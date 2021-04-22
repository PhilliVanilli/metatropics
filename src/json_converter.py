import csv
import pathlib

__author__ = 'Philippe Selhorst'


def json_converter(inputloc):
    """
    finds the sample_names.csv in the input location and converts it to rampart's json format into run_configuration.json
    """

    csvfile = pathlib.Path(inputloc, 'sample_names.csv')
    outputfile = pathlib.Path(inputloc, 'run_configuration.json')

    output = open(outputfile, 'w')
    output.write('{\n  "basecalledPath": "fastq/pass",\n  "samples": [\n')

    count = -1
    with open(csvfile, 'r') as handle:
        csv_reader = csv.reader(handle, dialect="excel")
        for line in csv_reader:
            count += 1

    with open(csvfile, 'r') as handle:
        csv_reader = csv.reader(handle, dialect="excel")
        for line_num, line in enumerate(csv_reader):
            if line_num != 0:
                barcode = line[0]
                barcode_number = barcode[-2:]
                new_barcode = 'NB' + barcode_number
                sample_name = line[2]
                json_block1 = '    {\n' + f'      "name": "{sample_name}",\n' + f'      "barcodes": [ "{new_barcode}" ]\n' + '    },\n'
                json_block2 = '    {\n' + f'      "name": "{sample_name}",\n' + f'      "barcodes": [ "{new_barcode}" ]\n' + '    }\n'
                if int(line_num) < int(count):
                    output.write(json_block1)
                else:
                    output.write(json_block2)
        output.write('  ]\n}\n')
        output.close()

