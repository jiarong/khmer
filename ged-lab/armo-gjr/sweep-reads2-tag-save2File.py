#! /usr/bin/env python
"""
Use a set of query reads to sweep out overlapping reads from another file.

% python scripts/sweep-reads2.py <query reads> <search reads> tag

Results end up in <search reads>.tag.sweep2, <search reads>.tag.nonswept.

Use '-h' for parameter help.
"""

import sys
import khmer
import os.path
import screed
from khmer.hashbits_args import build_construct_args, DEFAULT_MIN_HASHSIZE

def main():
    parser = build_construct_args()
    parser.add_argument('input_filename')
    parser.add_argument('read_filename')
    parser.add_argument('tag')

    args = parser.parse_args()

    if not args.quiet:
        if args.min_hashsize == DEFAULT_MIN_HASHSIZE:
            print >>sys.stderr, "** WARNING: hashsize is default!  " \
                "You absodefly want to increase this!\n** " \
                "Please read the docs!"

        print >>sys.stderr, '\nPARAMETERS:'
        print >>sys.stderr, ' - kmer size =    %d \t\t(-k)' % args.ksize
        print >>sys.stderr, ' - n hashes =     %d \t\t(-N)' % args.n_hashes
        print >>sys.stderr, ' - min hashsize = %-5.2g \t(-x)' % \
            args.min_hashsize
        print >>sys.stderr, ''
        print >>sys.stderr, 'Estimated memory usage is %.2g bytes ' \
            '(n_hashes x min_hashsize / 8)' % (
            args.n_hashes * args.min_hashsize / 8.)
        print >>sys.stderr, '-' * 8

    K = args.ksize
    HT_SIZE = args.min_hashsize
    N_HT = args.n_hashes

    inp = args.input_filename
    readsfile = args.read_filename
    tag = args.tag

    outfile = os.path.basename(readsfile) + '.' + tag + '.sweep2'
    nonswept = os.path.basename(readsfile) + '.' + tag + '.nonswept'
    outfp = open(outfile, 'wb')
    outfp2 = open(nonswept, 'wb')

    # create a hashbits data structure
    ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

    # load contigs, connect into N partitions
    print 'loading input reads from', inp
    ht.consume_fasta(inp)

    # Change 0.2 only if you really grok it.  HINT: You don't.
    fp_rate = khmer.calc_expected_collisions(ht)
    print >> sys.stderr, 'fp rate estimated to be %1.3f' % fp_rate

    if fp_rate > 0.20:
        print >>sys.stderr, "**"
        print >>sys.stderr, "** ERROR: the counting hash is too small for"
        print >>sys.stderr, "** this data set.  Increase hashsize/num ht."
        print >>sys.stderr, "**"
        print >>sys.stderr, "** Do not use these results!!"
        sys.exit(-1)

    print >> sys.stderr, 'starting sweep.'

    n = 0
    m = 0
    for record in screed.open(readsfile):
        if len(record.sequence) < K:
            continue

        if n % 1000000 == 0:
            print >> sys.stderr, '... %d %d' %(n, m)

        count = ht.get_median_count(record.sequence)[0]
        if count:
            m += 1
            outfp.write('>%s\n%s\n' % (record.name, record.sequence))
        else:
            outfp2.write('>%s\n%s\n' % (record.name, record.sequence))
        n += 1

    print >> sys.stderr, '%d of %d are swept' %(m, n)

if __name__ == '__main__':
    main()
