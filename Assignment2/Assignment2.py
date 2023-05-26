#! usr/bin/env python3

"""
Assignment2: program that receives fastq files and processes these
through a given server and multiprocessing returning the average
phredscores from the files.
"""

import multiprocessing as mp
from multiprocessing.managers import BaseManager, SyncManager
import os, sys, time, queue
import argparse as ap
import csv
import subprocess

POISONPILL = "MEMENTOMORI"
ERROR = "DOH"
IP = ''
PORTNUM = 1444
AUTHKEY = b'whathasitgotinitspocketsesss?'


def make_server_manager(port, authkey):
    """ Create a manager for the server, listening on the given port.
        Return a manager object with get_job_q and get_result_q methods.
    """
    job_q = queue.Queue()
    result_q = queue.Queue()

    # This is based on the examples in the official docs of multiprocessing.
    # get_{job|result}_q return synchronized proxies for the actual Queue
    # objects.
    class QueueManager(BaseManager):
        pass

    QueueManager.register('get_job_q', callable=lambda: job_q)
    QueueManager.register('get_result_q', callable=lambda: result_q)

    manager = QueueManager(address=('', port), authkey=authkey)
    manager.start()
    print('Server started at port %s' % port)
    return manager


def runserver(fn, data, sizes, output):
    # Start a shared manager server and access its queues
    manager = make_server_manager(PORTNUM, b'whathasitgotinitspocketsesss?')
    shared_job_q = manager.get_job_q()
    shared_result_q = manager.get_result_q()

    if not data:
        print("Gimme something to do here!")
        return

    print("Sending data!")
    for d in data:
        shared_job_q.put({'fn': fn, 'arg': d})

    time.sleep(2)

    results = []
    while True:
        try:
            result = shared_result_q.get_nowait()
            results.append(result)
            print("Got result!", result)
            if len(results) == len(data):
                print("Got all results!")
                break
        except queue.Empty:
            time.sleep(1)
            continue
    # Tell the client process no more data will be forthcoming
    print("Time to kill some peons!")
    shared_job_q.put(POISONPILL)
    # Sleep a bit before shutting down the server - to give clients time to
    # realize the job queue is empty and exit in an orderly way.
    time.sleep(5)
    print("Aaaaaand we're done for the server!")
    manager.shutdown()
    average_phredscores = calculate_average_phredscores(results, sizes)
    for file, scores in average_phredscores:
        if len(average_phredscores) > 1:
            if output is None:
                print(file)
                csvfile = None
            else:
                csvfile = f'{file}.{output}'
        else:
            csvfile = output
        create_output(scores, csvfile)



def make_client_manager(ip, port, authkey):
    """ Create a manager for a client. This manager connects to a server on the
        given address and exposes the get_job_q and get_result_q methods for
        accessing the shared queues from the server.
        Return a manager object.
    """
    class ServerQueueManager(BaseManager):
        pass

    ServerQueueManager.register('get_job_q')
    ServerQueueManager.register('get_result_q')

    manager = ServerQueueManager(address=(ip, port), authkey=authkey)
    manager.connect()

    print('Client connected to %s:%s' % (ip, port))
    return manager


def runclient(num_processes):
    manager = make_client_manager(IP, PORTNUM, AUTHKEY)
    job_q = manager.get_job_q()
    result_q = manager.get_result_q()
    run_workers(job_q, result_q, num_processes)


def run_workers(job_q, result_q, num_processes):
    processes = []
    for p in range(num_processes):
        temP = mp.Process(target=peon, args=(job_q, result_q))
        processes.append(temP)
        temP.start()
    print("Started %s workers!" % len(processes))
    for temP in processes:
        temP.join()


def peon(job_q, result_q):
    my_name = mp.current_process().name
    while True:
        try:
            job = job_q.get_nowait()
            if job == POISONPILL:
                job_q.put(POISONPILL)
                print("Aaaaaaargh", my_name)
                return
            else:
                try:
                    result = job['fn'](job['arg'])
                    print("Peon %s Workwork on %s!" % (my_name, job['arg']))
                    result_q.put({'job': job, 'result': result})
                except NameError:
                    print("Can't find yer fun Bob!")
                    result_q.put({'job': job, 'result': ERROR})

        except queue.Empty:
            print("sleepytime for", my_name)
            time.sleep(1)


def create_file_object(file, chunks_count):
    """
    Gets a fastq file and splits it in n amount of chunks, creating
    a start and end per chunk. Further it stores the size of the
    given file.
    :param file: fastq file
    :param chunks_count: Amount of chunks to split the file into
    :return: List object with chunks, file name and file size
    """
    quality_lines = subprocess.check_output(['wc', '-l', file])
    quality_lines = int(quality_lines.split()[0]) / 4
    chunks = []
    for i in range(chunks_count):
        start = round(i * quality_lines / chunks_count)
        end = round((i + 1) * quality_lines / chunks_count)
        chunks.append([file, start, end])
    return [chunks, file, quality_lines]


def read_fastq_chunk(chunk_object):
    """
    Reads and processes a chunk of a given fastq file and returns
    the quality scores of this chunk linked to its original file
    :param chunk_object: List containing a fastq file and a given
    start and end point to process from this file
    :return: Dictionary that links scores to its file
    """
    fastq_file = chunk_object[0]
    start = chunk_object[1]
    end = chunk_object[2]
    count = 0
    result = {}

    with open(fastq_file, encoding='UTF-8') as fastq:
        while count < start:
            fastq.readline()
            fastq.readline()
            fastq.readline()
            fastq.readline()
            count += 1
        while count < end:
            fastq.readline()
            fastq.readline()
            fastq.readline()
            quality = fastq.readline().rstrip()
            count += 1

            if not quality:
                # we reached the end of the file
                break

            for j, c in enumerate(quality):
                try:
                    result[fastq_file][j] += ord(c) - 33
                except KeyError:
                    result[fastq_file] = [ord(c) - 33]
                except IndexError:
                    result[fastq_file].append(ord(c) - 33)

    return result


def calculate_average_phredscores(results, num_reads):
    """
    Calculates average phredscores from the dictionarys received,
    adding scores to the linked file key
    :param results: Dictionary containing the input files as keys and
    a list with quality scores as values
    :param num_reads: Dictionary containing the input file as keys
    and the length of the files as values
    :return: Dictionary containing the average phredscores as values
    with the corresponding files as keys
    """
    phredscores = {}
    average_phredscores = {}
    for result in results:
        for file, scores in result.items():
            for i, score in enumerate(scores):
                try:
                    phredscores[file][i] += score
                except KeyError:
                    phredscores[file] = [score]
                except IndexError:
                    phredscores[file].append(score)
    for file, scores in phredscores.items():
        for score in scores:
            try:
                average_phredscores[file].append(score / num_reads[file])
            except KeyError:
                average_phredscores[file] = [score / num_reads[file]]

    return average_phredscores


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
    argparser = ap.ArgumentParser(
        description="Script voor Opdracht 2 van Big Data Computing;  Calculate PHRED scores over the network.")
    mode = argparser.add_mutually_exclusive_group(required=True)
    mode.add_argument("-s", action="store_true", help="Run the program in Server mode; see extra options needed below")
    mode.add_argument("-c", action="store_true", help="Run the program in Client mode; see extra options needed below")
    server_args = argparser.add_argument_group(title="Arguments when run in server mode")
    server_args.add_argument("-o", action="store", dest="csvfile", type=ap.FileType('w', encoding='UTF-8'),
                             required=False,
                             help="CSV file om de output in op te slaan. Default is output naar terminal STDOUT")
    server_args.add_argument("fastq_files", action="store", type=ap.FileType('r'), nargs='*',
                             help="Minstens 1 Illumina Fastq Format file om te verwerken")
    server_args.add_argument("--chunks", action="store", type=int, required=True)

    client_args = argparser.add_argument_group(title="Arguments when run in client mode")
    client_args.add_argument("-n", action="store",
                             dest="n", required=False, type=int,
                             help="Aantal cores om te gebruiken per host.")
    client_args.add_argument("--host", action="store", type=str, help="The hostname where the Server is listening")
    client_args.add_argument("--port", action="store", type=int, help="The port on which the Server is listening")

    args = argparser.parse_args()

    if args.s:
        if args.host:
            IP = args.host
        if args.port:
            PORTNUM = args.port
        if len(args.fastq_files) == 0:
            sys.exit("No input files received")
        file_objects = []
        jobs = []
        sizes = {}
        for file in args.fastq_files:
            file_objects.append(create_file_object(file, args.chunks))
        for obj in file_objects:
            for job in obj[0]:
                jobs.append(job)
            sizes[obj[1]] = obj[2]
        server = mp.Process(target=runserver, args=(read_fastq_chunk, jobs, sizes, args.csvfile))
        server.start()
        time.sleep(1)
        server.join()

    elif args.c:
        IP = args.host
        PORTNUM = args.port
        client = mp.Process(target=runclient, args=(args.n,))
        client.start()
        client.join()
    return 0


if __name__ == "__main__":
    sys.exit(main())
