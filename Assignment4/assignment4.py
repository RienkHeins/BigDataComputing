#! usr/bin/env python3

"""
Script that checks for the validity of a given fastq file and returns
the results to a csv file.
"""

import sys
import csv


def validate_file(fastq_file):
    """
    Gets a fasq file and checks for the validity by checking the line
    lengths and the starting characters of lines.
    :param fastq_file: fastq file to validate
    :return: validation information of input file
    """
    # open fastq file
    with open(fastq_file, encoding='UTF-8') as fastq:
        # create variables and parameters
        valid = True
        min_length = None
        max_length = 0
        counter = 0
        nuc_line_counter = 0
        total_length = 0
        # loop through file
        while True:
            # for each line check if the line is empty to see if end is reached
            # or if the file misses a line
            header = fastq.readline().rstrip()
            if len(header) == 0:
                break
            nucleotides = fastq.readline().rstrip()
            if len(nucleotides) == 0:
                counter += 1
                break
            strand = fastq.readline().rstrip()
            if len(strand) == 0:
                counter += 2
                break
            qual = fastq.readline().rstrip()
            if len(qual) == 0:
                counter += 3
                break
            # add the minimum length
            if min_length is None:
                min_length = len(nucleotides)
            # if line is shorter, set new minimum length
            elif min_length > len(nucleotides):
                min_length = len(nucleotides)
            # if line is longer, set new maximum length
            if max_length < len(nucleotides):
                max_length = len(nucleotides)
            # collect total length to calculate average
            total_length += len(nucleotides)
            nuc_line_counter += 1
            counter += 4
            # check if header, nucleotides and quality lines are correct
            if valid:
                if not header.startswith("@"):
                    valid = False
                if len(nucleotides) != len(qual):
                    valid = False
        # check if the file doesn't miss any lines
        if counter % 4 != 0:
            valid = False
        return [fastq_file, valid, min_length, max_length, total_length/nuc_line_counter]


def create_output(row):
    """
    Creates csv output from the input file and writes it to the
    command line output
    :param row: list object with the filename, minimal line length,
    maximal line length and average line length.
    """
    # generate output
    header = ["Filenaam,Valide", "Min_length", "Max_length", "Average_length"]
    csv.writer(sys.stdout, delimiter=",").writerow(header)
    csv.writer(sys.stdout, delimiter=",").writerow(row)


def main(file):
    """
    main file to generate output
    :param file: fastq file
    """
    create_output(validate_file(file))


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
