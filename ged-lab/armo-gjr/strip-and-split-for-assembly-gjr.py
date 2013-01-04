#! /usr/bin/env python
# this use big mem if most of seqs in a big file can pair up.
# by gjr; Nov 30, 12
import sys
import screed
from screed.fasta import fasta_iter
import random

fwSE = open('%s.se' %(sys.argv[1]), 'wb')
fwPE = open('%s.pe' %(sys.argv[1]), 'wb')
d = {}
for n, record in enumerate(screed.open(sys.argv[1])):
    name = record.name
    if name.endswith('/1') or name.endswith('/2'):
        key = name.rstrip('/1').rstrip('/2')
        x = d.get(key, [])
        x.append(name)
        d[key] = x
    else:
        seq = record.sequence
        print >> fwSE, '>%s\n%s' %(name, seq)

cnt = 0
dSeq = {}
for record in screed.open(sys.argv[1]):
    name = record.name
    if name.endswith('/1') or name.endswith('/2'):
        key = name.rstrip('/1').rstrip('/2')
        seq = record.sequence
        if len(d[key]) == 1:
            x = d.pop(key)
            print >> fwSE, '>%s\n%s' %(name, seq)
        else:
            dSeq[name] = seq   ### may take big mem if most seqs are pairEnded
            cnt += 1
            

assert 2*len(d) == len(dSeq), 'should equals to number of pairs'
print '%d pairs are found in %s' %(cnt/2, sys.argv[1])
for key in d:
    name1, name2 = sorted(d[key])
    print >> fwPE, '>%s\n%s' %(name1, dSeq[name1])
    print >> fwPE, '>%s\n%s' %(name2, dSeq[name2])

