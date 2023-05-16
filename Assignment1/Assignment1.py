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
    splitted = []
    size = len(data)
    for i in range(chunks):
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

    with open(fastq_file, encoding='UTF-8') as fastq:
        while quality:
            _ = fastq.readline()
            _ = fastq.readline()
            _ = fastq.readline()
            quality = fastq.readline().rstrip()

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
    for qual in quals:
        for i, char in enumerate(qual):
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
    if csvfile is None:
        csv_writer = csv.writer(sys.stdout, delimiter=',')
        for i, score in enumerate(average_phredscores):
            csv_writer.writerow([i, score])

    else:
        with open(csvfile, 'w', encoding='UTF-8', newline='') as myfastq:
            csv_writer = csv.writer(myfastq, delimiter=',')
            for i, score in enumerate(average_phredscores):
                csv_writer.writerow([i, score])


def main():
    argparser = ap.ArgumentParser(description="Script for assignment 1 of Big Data Computing")
    argparser.add_argument("-n", action="store",
                           dest="n", required=True, type=int,
                           help="Amount of cores to be used")
    argparser.add_argument("-o", action="store", dest="csvfile", required=False,
                           help="CSV file to save the output. Default is output to terminal STDOUT")
    argparser.add_argument("fastq_files", action="store",
                           nargs='+', help="At least 1 ILLUMINA fastq file to process")
    args = argparser.parse_args()

    for file in args.fastq_files:
        qualities = read_fastq_file(file)
        qualities_chunked = chunks(qualities, 4)

        with mp.Pool(args.n) as pool:
            phredscores = pool.map(calculate_quals, qualities_chunked)
        phredscores_avg = [sum(i) / len(qualities) for i in zip(*phredscores)]
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
