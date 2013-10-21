#! /usr/bin/env python

import khmer
import sys, threading, time, argparse, cPickle
from khmer.threading_args import add_threading_args
 
"""
Count the median/avg k-mer abundance for each sequence in the input file,
based on the k-mer counts in the given counting hash.  Can be used to
estimate expression levels (mRNAseq) or coverage (genomic/metagenomic).

% scripts/count-median.py <htname> <input seqs> <output counts>

Use '-h' for parameter help.

The output is pickled dict contains sequence id, median.

NOTE: All 'N's in the input sequences are converted to 'G's.
"""

#def main():
start = time.time()
parser = argparse.ArgumentParser(
    description='Count k-mers summary stats for sequences')

add_threading_args(parser)

parser.add_argument('htfile')
parser.add_argument('input')
parser.add_argument('output')

args = parser.parse_args()

htfile = args.htfile
input_filename = args.input
output_filename = args.output

n_threads = int(args.n_threads)

config = khmer.get_config()
bufsz = config.get_reads_input_buffer_size()
#default_threads = config.get_number_of_threads()
#print '>>>>> bufsz: %d; default_threads: %d' %(bufsz, default_threads)

config.set_number_of_threads(n_threads)
new_bufsz = n_threads * bufsz
config.set_reads_input_buffer_size(new_bufsz)
rparser = khmer.ReadParser(input_filename, n_threads)
print >> sys.stderr, '### buffer size: %d; threads: %d' %(new_bufsz,\
                                                            n_threads)

print >> sys.stderr, 'loading counting hash from %s' %htfile
ht = khmer.load_counting_hash(htfile)
end1 = time.time()
print >> sys.stderr, 'loading took %d sec' %(end1-start)
K = ht.ksize()
COUNTING_DICT = {}

for i in xrange(n_threads):
    # initial dict of dict
    COUNTING_DICT[i] = {}


###
def count_median(rparser, tnum):
    for n, record in enumerate(rparser):
        seq = record.sequence.upper()
        name = record.name
        if 'N' in seq:
            seq = seq.replace('N', 'G')


        if K <= len(seq):
            a, b, c = ht.get_median_count(seq)
            COUNTING_DICT[i][name] = a

        if n%1e6 == 0:
            print '%d process by thread %d' %(n, tnum)
###
print >> sys.stderr, 'writing to %s' %output_filename
output = open(output_filename, 'wb')


threads = []

for tnum in xrange(n_threads):
    print >> sys.stderr, 'start counting with %d threads' %tnum
    t = threading.Thread(
            target=count_median,
            args=(rparser,tnum,)
    )
    threads.append(t)
    t.start()

for t in threads:
    t.join()

end2 = time.time()

MERGED_DICT = {}
for i in xrange(n_threads):
    di = COUNTING_DICT[i]
    MERGED_DICT.update(di)

print >> sys.stderr, 'parsing took %d sec' %(end2-end1)
print >> sys.stderr, 'finshed..'
print >> sys.stderr, '%d reads processed in total' %len(MERGED_DICT)
cPickle.dump(MERGED_DICT, open('%s' %output_filename, 'wb'))

#if __name__ == '__main__':
#    main()
