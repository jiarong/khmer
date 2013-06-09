#! /usr/bin/python
# convert fastq to fasta
# by gjr; Oct 4, 11


import sys
sys.path.insert(0, '/mnt/home/guojiaro/Documents/lib/screed')
import screed

f = sys.argv[1]
fw = open(sys.argv[2], 'w')

for n, record in enumerate(screed.open(f)):
    name = record['name']
    seq = record['sequence']
    print >> fw, '>%s\n%s' %(name, seq)

print (n+1), 'fasta seqs written'
fw.close()
