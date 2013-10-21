#! /usr/bin/env python
import screed
import sys

CUTOFF = 300


n = 0
for filename in sys.argv[1:]:
    print >> sys.stderr, 'processing %s' %(filename)
    for record in screed.open(filename):
        if len(record.sequence) >= CUTOFF:
            n += 1
            print '>%s\n%s' % (n, record.sequence)
