#! usr/bin/env python3

"""
Multiprocessing program that reads and processes fastq files and
returns a csv file with the average phred scores of these files.
"""

import sys
import argparse as ap
import multiprocessing as mp
import csv


def chunks(data, chunks):
    """
    Divides the data into multiple chunks for multiprocessing
    :param data: Quality scores from fastq files
    :param chunks: Amount of chunks to divide the data into
    :return: Splitted data
    """
    # create list for splitted data and get the size of the file
    splitted = []
    size = len(data)
    # create n amount of chunks
    for i in range(chunks):
        # split the data based on chunks and file size
        start = int(i*size/chunks)
        end = int((i+1)*size/chunks)
        splitted.append(data[start:end])
    return splitted


def read_fastq_file(fastq_file):
    """
    Reads a fastq file and returns the quality scores from this file
    :param fastq_file: fastq file
    :return: quality scores from the given file
    """
    quality_scores = []
    quality = True
    # open file and loop through
    with open(fastq_file, encoding='UTF-8') as fastq:
        # check if quality contains characters, if not end of the file is reached
        while quality:
            # skipp files without quality information
            fastq.readline()
            fastq.readline()
            fastq.readline()
            # strip line
            quality = fastq.readline().rstrip()
            # append if it contains info
            if quality:
                quality_scores.append(quality)
    return quality_scores


def calculate_quals(quals):
    """
    Calculates quality scores
    :param quals: list of fastq quality score lines
    :return: quality scores
    """
    results = []
    # loop through quality lines
    for qual in quals:
        for i, char in enumerate(qual):
            # calculate scores
            try:
                results[i] += ord(char) - 33
            except IndexError:
                results.append(ord(char) - 33)
    return results


def create_output(average_phredscores, csvfile):
    """
    Generates csv output files
    :param average_phredscores: average phredscores calculated from
    the fastq files
    :param csvfile: csv output file name
    """
    # check if a filename is given
    if csvfile is None:
        # write results
        csv_writer = csv.writer(sys.stdout, delimiter=',')
        for i, score in enumerate(average_phredscores):
            csv_writer.writerow([i, score])
    else:
        # write results
        with open(csvfile, 'w', encoding='UTF-8', newline='') as myfastq:
            csv_writer = csv.writer(myfastq, delimiter=',')
            for i, score in enumerate(average_phredscores):
                csv_writer.writerow([i, score])


def main():
    # create argparser
    argparser = ap.ArgumentParser(description="Script for assignment 1 of Big Data Computing")
    argparser.add_argument("-n", action="store",
                           dest="n", required=True, type=int,
                           help="Amount of cores to be used")
    argparser.add_argument("-o", action="store", dest="csvfile", required=False,
                           help="CSV file to save the output. Default is output to terminal STDOUT")
    argparser.add_argument("fastq_files", action="store",
                           nargs='+', help="At least 1 ILLUMINA fastq file to process")
    args = argparser.parse_args()
    # loop through files
    for file in args.fastq_files:
        qualities = read_fastq_file(file)
        qualities_chunked = chunks(qualities, 4)
        # create multiprocessing pools
        with mp.Pool(args.n) as pool:
            phredscores = pool.map(calculate_quals, qualities_chunked)
        # calculate average phredscores
        phredscores_avg = [sum(i) / len(qualities) for i in zip(*phredscores)]
        # write output
        if len(args.fastq_files) > 1:
            if args.csvfile is None:
                print(file)
                csvfile = None
            else:
                csvfile = f'{file}.{args.csvfile}'
        else:
            csvfile = args.csvfile
        create_output(phredscores_avg, csvfile)


if __name__ == "__main__":
    sys.exit(main())
