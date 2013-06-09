#! /usr/bin/env python
import screed
import sys
import os.path

fp1 = open(os.path.basename(sys.argv[1]) + '.1', 'wb')
fp2 = open(os.path.basename(sys.argv[1]) + '.2', 'wb')

n1 = 0
n2 = 0
for n, record in enumerate(screed.open(sys.argv[1])):
    if n % 1000000 == 0:
        print >>sys.stderr, '...', n

    anno = record.annotations
    if anno.startswith('1:'):
    #if name.endswith('.1'):
        #print >>fp1, '>%s\n%s' % (record.name, record.sequence,)
        print >>fp1, '@%s/1\n%s\n+\n%s' % (record.name, record.sequence, record.accuracy)
        n1 += 1
    elif anno.startswith('2:'):
    #elif name.endswith('.2'):
        print >>fp2, '@%s/2\n%s\n+\n%s' % (record.name, record.sequence, record.accuracy)
        n2 += 1

print >>sys.stderr, "DONE; split %d sequences (%d left, %d right)" % \
    (n + 1, n1, n2)
