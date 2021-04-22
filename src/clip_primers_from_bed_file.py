import argparse
from json import dump
import pathlib
from copy import copy
import csv
from operator import itemgetter
import pysam
import collections

"""
Written by Nick Loman as part of the ZiBRA pipeline (zibraproject.org)
edited by Colin Anthony
"""


def check_still_matching_bases(s):
    for flag, length in s.cigartuples:
        if flag == 0:
            return True
    return False


def read_bed_file(primerset_bed):

    out_bedfile = []
    with open(primerset_bed) as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel-tab')
        for row in reader:
            primer_name = row['Primer_ID']
            orientation = primer_name.split("_")
            if "LEFT" in orientation:
                row['direction'] = "+"
            elif "RIGHT" in orientation:
                row['direction'] = "-"
            else:
                print("could not parse primer name\nexpected either LEFT or RIGHT in name, using '_' as a separator")
                raise ValueError

            row['end'] = int(row['end'])
            row['start'] = int(row['start'])
            out_bedfile.append(row)

    return out_bedfile


def trim(cigar, s, start_pos, end):
    if not end:
        pos = s.pos
    else:
        pos = s.reference_end

    eaten = 0
    while 1:
        # chomp stuff off until we reach pos
        if end:
            try:
                flag, length = cigar.pop()
            except IndexError:
                total_len = sum([x[1] for x in s.cigartuples])
                soft_clip = sum([x[1] for x in s.cigartuples if x[0] == 4])

                print(f"Nothing left after cipping\n"
                      f"most likely wrong primer scheme used\n"
                      f"sequence name = {s.query_name}\n"
                      f"cigar code={s.cigarstring}\n"
                      f"total length = {total_len}\n"
                      f"total_soft_clipped_by_bwa={soft_clip}\n"
                      f"total clipped = {eaten}\n")
                break
        else:
            try:
                flag, length = cigar.pop(0)
            except IndexError:
                total_len = sum([x[1] for x in s.cigartuples])
                soft_clip = sum([x[1] for x in s.cigartuples if x[0] == 4])

                print(f"Nothing left after cipping\n"
                      f"most likely wrong primer scheme used\n"
                      f"sequence name = {s.query_name}\n"
                      f"cigar code={s.cigarstring}\n"
                      f"total length = {total_len}\n"
                      f"total_soft_clipped_by_bwa={soft_clip}\n"
                      f"total clipped = {eaten}\n")
                break

        if flag == 0:
            # match
            eaten += length
            if not end:
                pos += length
            else:
                pos -= length
        if flag == 1:
            # insertion to the ref
            eaten += length
        if flag == 2:
            # deletion to the ref
            if not end:
                pos += length
            else:
                pos -= length
            pass
        if flag == 4:
            # soft clip
            eaten += length
        if not end and pos >= start_pos and flag == 0:
            break
        if end and pos <= start_pos and flag == 0:
            break

    extra = abs(pos - start_pos)

    if extra:
        if flag == 0:
            if end:
                cigar.append((0, extra))
            else:
                cigar.insert(0, (0, extra))
            eaten -= extra

    if not end:
        s.pos = pos - extra
    if end:
        cigar.append((4, eaten))
    else:
        cigar.insert(0, (4, eaten))

    if cigar[0][1] <= 0 or cigar[-1][1] <= 0:
        print("negative length added to cigar, probable indel in primer region")
        # print("old", s.cigartuples)
        # print("new", cigar)
        raise
    s.cigartuples = cigar

    return cigar


def find_primer(bed, ref_pos, direction):
    closest = min([(abs(p['start'] - ref_pos), p['start'] - ref_pos, p) for p in bed if p['direction'] == direction],
                  key=itemgetter(0))
    return closest


def is_correctly_paired(p1, p2):
    name1 = p1[2]['Primer_ID']
    name2 = p2[2]['Primer_ID']

    name1 = name1.replace('_LEFT', '')
    name2 = name2.replace('_RIGHT', '')

    return name1 == name2


def main(infile, outfile, bedfile):

    bed = read_bed_file(bedfile)
    sam_infile = pysam.AlignmentFile(infile, "r")
    outfile = pathlib.Path(outfile).absolute()
    outfile_trimmed = pysam.AlignmentFile(str(outfile), "wh", template=sam_infile)
    suppl_out = str(outfile) + "_excluded_sequences_as_Supplementary.sam"
    marked_supplementary = pysam.AlignmentFile(suppl_out, "wh", template=sam_infile)
    primer_mismatch_file = str(outfile) + "_excluded_as_primer_mismatched.sam"
    marked_primer_missmatch = pysam.AlignmentFile(primer_mismatch_file, "wh", template=sam_infile)

    runfolder = pathlib.Path(outfile.parent.parent.parent)
    runname = runfolder.parts[-1]
    mappingfile = pathlib.Path(runfolder, runname + '_mapping.csv')

    # set counters
    total = 0
    unmapped = 0
    missmatched = 0
    suppl = 0
    good = 0
    bad = 0

    read_prime_pair_lookup_dict = collections.defaultdict(list)

    # make sure that all primer pairs are represented in the dict
    for bed_index in range(0, len(bed), 2):
        if bed[bed_index]["number"] != bed[bed_index + 1]["number"]:
            print("primers are not a pair")
            raise ValueError

        if bed[bed_index]["direction"] == "+":
            primer_start = bed[bed_index]["end"]
        else:
            raise ValueError
        if bed[bed_index + 1]["direction"] == "-":
            primer_end = bed[bed_index + 1]["start"]
        else:
            raise ValueError
        pair_key = f"{primer_start}_{primer_end}"
        read_prime_pair_lookup_dict[pair_key] = []

    for s in sam_infile:
        total += 1
        cigar = copy(s.cigartuples)

        # logic - if alignment start site is _before_ but within X bases of  a primer site, trim it off
        if s.is_unmapped:
            # print(f"{s.query_name} skipped as unmapped")
            unmapped += 1
            continue

        # logic to remove part of read that were mapped to non-contiguous region (ie supplementary read)
        if s.is_supplementary:
            marked_supplementary.write(s)
            suppl += 1
            continue

        p1 = find_primer(bed, s.reference_start, '+')
        p2 = find_primer(bed, s.reference_end, '-')

        if not is_correctly_paired(p1, p2):
            # print("mismatched primer pair. primers matched:", p1[2]['Primer_ID'], p2[2]['Primer_ID'],
            #       "this is probably two amplicons ligated together")
            marked_primer_missmatch.write(s)
            missmatched += 1
            continue

        # create dict to write seq name and primer pair code to json

        read_prime_pair_lookup_dict[f"{p1[2]['end']}_{p2[2]['start']}"].append(s.query_name)
        # if the alignment starts before the end of the primer, trim to that position
        primer_position = p1[2]['end']
        pass_1 = False
        pass_2 = False
        try:
            if s.reference_start < primer_position:
                trim(cigar, s, primer_position, 0)
                pass_1 = True
        except Exception as e:
            print("problem clipping primers, most likely due to indels in the primer region\n skipping read\n")
            bad += 1
            continue

        primer_position = p2[2]['start']

        try:
            if s.reference_end > primer_position:
                trim(cigar, s, primer_position, 1)
                pass_2 = True
        except Exception as e:
            print("problem clipping primers, most likely due to indels in the primer region\n skipping read\n")
            bad += 1
            continue

        if pass_1 and pass_2:
            good += 1

        if not check_still_matching_bases(s):
            continue

        outfile_trimmed.write(s)

    # write seq name and primer pair code to json
    with open(read_prime_pair_lookup, 'w') as jd:
        dump(read_prime_pair_lookup_dict, jd)

    print(f"Total: {total}\nUnmapped: {unmapped} ({round(unmapped/total*100, 2)}%)\n"
          f"Supplementary: {suppl} ({round(suppl/total*100, 2)}%)\n"
          f"Mismatched primers: {missmatched} ({round(missmatched/total*100, 2)}%)\n"
          f"indel in primer sequences: {bad} ({round(bad/total*100, 2)}%)\n"
          f"Good sequences: {good} ({round(good/total*100, 2)}%)\n")

    with open(mappingfile,'a') as fh:
        fh.write(f"\n{outfile.stem}\n"
          f"Total,{total},100%\nUnmapped,{unmapped},{round(unmapped/total*100, 2)}%\n"
          f"Supplementary,{suppl},{round(suppl/total*100, 2)}%\n"
          f"Mismatched primers,{missmatched},{round(missmatched/total*100, 2)}%\n"
          f"indel in primer sequences,{bad},{round(bad/total*100, 2)}%\n"
          f"Good sequences,{good},{round(good/total*100, 2)}%\n")




    print("\nFinished soft clipping bam file\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Trim alignments from an amplicon scheme.')
    parser.add_argument('-in', '--infile', help='Input BAM filename and path')
    parser.add_argument('-o', '--outfile', help='output SAM filename and path')
    parser.add_argument('-b', '--bedfile', help='BED file containing the amplicon scheme')

    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile
    bedfile = args.bedfile

    main(infile, outfile, bedfile)
