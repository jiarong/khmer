#! /usr/bin/env python
import screed
import sys
import os.path

#usage: python <thisfile> <inputNewFastq> <outputOldfastq>

try:
    input = sys.argv[1]
    output = sys.argv[2]
except IndexError:
    print 'usage: python <thisfile> <inputNewFastq> <outputOldfastq>'
    sys.exit(1)

fw = open(output, 'wb')
for n, record in enumerate(screed.open(input)):
    if n % 1000000 == 0:
        print >>sys.stderr, '...', n

    anno = record.annotations
    if anno.startswith('1:'):
        print >>fw, '@%s/1\n%s\n+\n%s' % (record.name, record.sequence, record.accuracy)
    elif anno.startswith('2:'):
        print >>fw, '@%s/2\n%s\n+\n%s' % (record.name, record.sequence, record.accuracy)

print >>sys.stderr, "DONE;  processed %d sequences" % \
    (n + 1)
