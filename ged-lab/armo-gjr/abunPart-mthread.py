#! /usr/bin/env python

import khmer
import sys, threading, time, argparse, cPickle, math, os
import array
import subprocess
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
parser.add_argument('label')
parser.add_argument('upcutoff')
parser.add_argument('lowcutoff')

args = parser.parse_args()

htfile = args.htfile
input_filename = args.input
label = args.label

NUM = int(args.upcutoff)
NUM = math.log(NUM, 2)
NUM = int(math.floor(NUM))
NUM = 2**NUM

lowcutoff = int(args.lowcutoff)

n_threads = int(args.n_threads)
print >> sys.stderr, 'threads used: %d' %n_threads

config = khmer.get_config()
bufsz = config.get_reads_input_buffer_size()
config.set_reads_input_buffer_size(n_threads * 1 * 1024 * 1024)
rparser = khmer.ReadParser(input_filename, n_threads)

print >> sys.stderr, 'loading counting hash from %s' %htfile
ht = khmer.load_counting_hash(htfile)
K = ht.ksize()
MAXCOUNT = 2**16 - 1
COUNTING_LIS = [0]*(MAXCOUNT+1)
DICT_CNT_ARRAY = {}
DICT_MEAN1_CNT = {}
DICT_MED1_FILT_CNT = {}
DICT_MEAN1_FILE = {}
DICT_MEAN1_FILE_NAMES = {}
DICT_MED1_FILE = {}
DICT_MED1_FILE_NAMES = {}
dF = {}
dF_NAMES = {}

for tnum in xrange(n_threads):
    DICT_CNT_ARRAY[tnum] = array.array('I', COUNTING_LIS)
    DICT_MEAN1_CNT[tnum] = 0
    DICT_MED1_FILT_CNT[tnum] = 0
    fname = '%s.mean1.fasta.thread%d' %(label, tnum)
    DICT_MEAN1_FILE_NAMES[tnum] = fname

    DICT_MEAN1_FILE[tnum] = open(fname, 'wb')
    fname = '%s.%d.filtered.fasta.thread%d' %(label,lowcutoff,tnum)
    DICT_MED1_FILE_NAMES[tnum] = fname
    DICT_MED1_FILE[tnum] = open(fname, 'wb')
    
    i = 2
    dF[tnum] = {}
    while (i <= NUM):
        if (i*4 <= lowcutoff):
            i *= 2
            continue
        low_bound = i
        if (i < lowcutoff):
            low_bound = lowcutoff+1
        fname = '%s.%dto%d.fasta.thread%d' %(label, low_bound, i*4, tnum)
        dF_NAMES.setdefault(i, [])
        dF_NAMES[i].append(fname)
        dF[tnum][i] = open(fname, 'wb')
        i *= 2
    fname = '%s.%dtoMAX.fasta.thread%d' %(label, i, tnum)
    dF_NAMES.setdefault(i, [])
    dF_NAMES[i].append(fname)
    dF[tnum][i] = open(fname, 'wb')

end1 = time.time()
print >> sys.stderr, 'loading took: %d sec' %(end1 - start)
###
def count_median(rparser, tnum):
    for n, record in enumerate(rparser):
        seq = record.sequence.upper()
        name = record.name
        if 'N' in seq:
            seq = seq.replace('N', 'G')

        if K <= len(seq):
            a, b, c = ht.get_median_count(seq)
            if a > MAXCOUNT:
                a = MAXCOUNT
            DICT_CNT_ARRAY[tnum][a] += 1
            if (a <= lowcutoff):
                if (b == 1):
                    print >> DICT_MEAN1_FILE[tnum], '>%s\n%s' %(record.name,\
                                                                  seq)
                    DICT_MEAN1_CNT[tnum] += 1
                else:
                    print >> DICT_MED1_FILE[tnum], '>%s\n%s' %(record.name,\
                                                                   seq)
                    DICT_MED1_FILT_CNT[tnum] += 1
            else:
                tempN = math.log(a, 2)
                tempN = int(math.floor(tempN))
                tempN = 2**tempN
                if tempN >= NUM*4:
                    print >> dF[tnum][NUM*2], '>%s\n%s' %(record.name, seq)
                else:
                    if tempN/2 != 1: 
                        print >> dF[tnum][tempN/2], \
                                    '>%s\n%s' %(record.name, seq)
                    print >> dF[tnum][tempN], '>%s\n%s' %(record.name, seq)

        if n%1e6 == 0:
            print >> sys.stderr, '%d process by thread %d' %(n, tnum)
###


threads = []

for tnum in xrange(n_threads):
    #print >> sys.stderr, 'start counting with %d threads' %tnum
    t = threading.Thread(
            target=count_median,
            args=(rparser,tnum,)
    )
    threads.append(t)
    t.start()

for t in threads:
    t.join()

end2 = time.time()
# wait for all threads to print to std
time.sleep(10)
print >> sys.stderr, 'counting took: %d sec' %(end2 - end1)


# merge array
merged_array = zip(*DICT_CNT_ARRAY.values())
for i in xrange(len(COUNTING_LIS)):
    COUNTING_LIS[i] = sum(merged_array[i])

# add up variables
mean1_cnt = sum(DICT_MEAN1_CNT.values())
med1_filt_cnt = sum(DICT_MED1_FILT_CNT.values())

print >> sys.stderr, 'finshed..'
print >> sys.stderr, '%d reads processed in total' %sum(COUNTING_LIS)
print >> sys.stderr, 'number of reads with mean kmer count 1 (singleton): %d'\
                          %mean1_cnt
print >> sys.stderr, 'number of reads with median kmer count <= %d' \
                          '(not including singleton): %d' \
                          %(lowcutoff, med1_filt_cnt)

cnt = len(COUNTING_LIS)-1
while (cnt >= 0):
    if (COUNTING_LIS[cnt] != 0):
        max_med_cnt = cnt
        max_med_cnt_abun = COUNTING_LIS[cnt]
        break
    cnt -= 1
    
print >> sys.stderr, 'max median count is %d with %d reads  ' \
                        %(max_med_cnt, max_med_cnt_abun)

# close files to flush
for tnum in xrange(n_threads):
    DICT_MEAN1_FILE[tnum].close()
    DICT_MED1_FILE[tnum].close()
    for fp in dF[tnum].values():
        fp.close()

# merge files
mean1_file_name_lis = DICT_MEAN1_FILE_NAMES.values()
mean1_file_name = mean1_file_name_lis[0].rsplit('.', 1)[0]
p = subprocess.Popen('cat %s > %s' % (' '.join(mean1_file_name_lis), 
                                                  mean1_file_name), 
                                      shell=True,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
(stdout, stderr) = p.communicate()
assert p.returncode == 0, (stdout, stderr)
#remove .thread files
for fname in mean1_file_name_lis:
    os.remove(fname)

# merge
med1_file_name_lis = DICT_MED1_FILE_NAMES.values()
med1_file_name = med1_file_name_lis[0].rsplit('.', 1)[0]
p = subprocess.Popen('cat %s > %s' % (' '.join(med1_file_name_lis),
                                                  med1_file_name),
                                      shell=True,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
(stdout, stderr) = p.communicate()
assert p.returncode == 0, (stdout, stderr)

#remove .thread files
for fname in med1_file_name_lis:
    os.remove(fname)


for flis in dF_NAMES.values():
    assert flis[0].rsplit('.',1)[0] == flis[-1].rsplit('.',1)[0]
    #merge
    file_name = flis[0].rsplit('.', 1)[0]
    p = subprocess.Popen('cat %s > %s' % (' '.join(flis),
                                                      file_name),
                                          shell=True,
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE)
    (stdout, stderr) = p.communicate()
    assert p.returncode == 0, (stdout, stderr)
    #remove .thread files
    for fname in flis:
        os.remove(fname)
    
tempN = math.log(max_med_cnt, 2)
tempN = int(math.floor(tempN))+1

NUM_log = int(math.log(NUM,2))

if tempN < NUM_log + 1:
    os.remove('%s.%dtoMAX.fasta' %(label, NUM*2))
    for i in range(tempN, NUM_log+1):
        ind = 2**i
        os.remove('%s.%dto%d.fasta' %(label, ind, ind*4))


output = open('%s.dist' %label, 'wb')
for i in xrange(len(COUNTING_LIS)):
    print >> output, '%d\t%d' %(i, COUNTING_LIS[i])


#if __name__ == '__main__':
#    main()
