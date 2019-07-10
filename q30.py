#!/usr/bin/env python

import os,sys
import re
import fastq
import time

def qual_stat(qstr):
    q20 = 0
    q30 = 0
    for q in qstr:
        qual = ord(q) - 33
        if qual >= 30:
            q30 += 1
            q20 += 1
        elif qual >= 20:
            q20 += 1
    return q20, q30

def gc_stat(seqs):
    total = 0
    for i in seqs:
        if i in ['G', 'C']:
            total += 1
    return total

def seq_len(seqs):
    return len(seqs)

def read_count(filename):
    if re.search(r'fa$', filename, re.I):
        reads_count = 0
        with open(filename) as f:
            for line in f:
                if line.startswith('>'):
                    reads_count += 1
        return reads_count
    else:
        reader = fastq.Reader(filename)
        reads_count = 0
        while True:
            read = reader.nextRead()
            if read == None:
                break
            reads_count += 1
        return reads_count
    # print(reads_count)


def stat(filename):
    reader = fastq.Reader(filename)
    total_count = 0
    q20_count = 0
    q30_count = 0
    gc_count = 0
    length_count = 0
    reads_count = 0
    while True:
        read = reader.nextRead()
        if read == None:
            break
        total_count += len(read[3])
        q20, q30 = qual_stat(read[3])
        q20_count += q20
        q30_count += q30
        gc = gc_stat(read[1])
        gc_count += gc
        length = seq_len(read[1])
        length_count += length
        reads_count += 1

    q20_percents = round(100 * float(q20_count)/float(total_count), 2)
    q30_percents = round(100 * float(q30_count)/float(total_count), 2)
    GC_percents = round(100 * float(gc_count)/float(total_count), 2)
    average_len = round(float(length_count)/float(reads_count))

    return total_count, average_len, q20_percents, q30_percents, GC_percents

    # print(reads_count)
    # print("total bases:", total_count)
    # print("q20 bases:", q20_count)
    # print("q30 bases:", q30_count)
    # print('GC bases:', gc_count)
    # print("q20 percents:", round(100 * float(q20_count)/float(total_count), 2))
    # print("q30 percents:", round(100 * float(q30_count)/float(total_count), 2))
    # print('GC percents:', round(100 * float(gc_count)/float(total_count), 2))
    # print('average length:', round(float(length_count)/float(reads_count)))

def main():
    if len(sys.argv) < 2:
        print("usage: python q30.py <fastq_file>")
        sys.exit(1)
    stat(sys.argv[1])
    read_count(sys.argv[1])

if __name__ == "__main__":
    time1 = time.time()
    main()
    time2 = time.time()
    print('Time used: ' + str(time2-time1))