import csv
import pathlib

__author__ = 'Philippe Selhorst'


def json_converter(inputloc, barcodes):
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
                if barcodes == 'CUST':
                    prefix = ('BC')
                    cb_dict = {'01':'04', '02':'06', '03':'09', '04':'12', '05':'18', '06':'19', '07':'20', '08':'21',
                              '09':'22', '10':'31', '11':'38', '12':'46', '13':'50', '14':'55', '15':'56', '16':'67',
                              '17':'72', '18':'75', '19':'78', '20':'79', '21':'80', '22':'81', '23':'86', '24':'88',
                              '25':'89', '26': '95', '27':'96'}
                    new_barcode = prefix + cb_dict[barcode_number]
                else:
                    prefix = ('NB')
                    new_barcode = prefix + barcode_number
                sample_name = line[2]
                json_block1 = '    {\n' + f'      "name": "{sample_name}",\n' + f'      "barcodes": [ "{new_barcode}" ]\n' + '    },\n'
                json_block2 = '    {\n' + f'      "name": "{sample_name}",\n' + f'      "barcodes": [ "{new_barcode}" ]\n' + '    }\n'
                if int(line_num) < int(count):
                    output.write(json_block1)
                else:
                    output.write(json_block2)
        output.write('  ]\n}\n')
        output.close()

