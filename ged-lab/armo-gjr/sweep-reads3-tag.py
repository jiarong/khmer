#! /usr/bin/env python
"""
Use a set of query reads to sweep out overlapping reads from multiple files.

% python scripts/sweep-reads3.py <query1> [ <query2> ... ] <search reads> tag

Results end up in <query?>.tag.sweep3.

Use '-h' for parameter help.
"""


import sys
import os.path
import screed
import khmer
from khmer.hashbits_args import build_construct_args, DEFAULT_MIN_HASHSIZE

def main():
    parser = build_construct_args()
    parser.add_argument('input_filenames', nargs='+')
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
            args.n_hashes * args.min_hashsize * len(args.input_filenames) / 8.)
        print >>sys.stderr, '-'*8

    K=args.ksize
    HT_SIZE=args.min_hashsize
    N_HT=args.n_hashes

    inputlist = args.input_filenames 
    readsfile = args.read_filename
    tag = args.tag

    logfile = open('%s.sweep3.log' %(tag), 'wb')
    query_list = []
    for n, inp_name in enumerate(inputlist):
        # create a hashbits data structure
        ht = khmer.new_hashbits(K, HT_SIZE, N_HT)

        outfile = os.path.basename(inp_name) + '.' + tag + '.sweep3'
        outfp = open(outfile, 'w')
        query_list.append((ht, outfp))

    new_lis = []
    cnt = 0
    for n, inp_name in enumerate(inputlist):
        ht = query_list[n][0]
        outfp = query_list[n][1]

        # load contigs, connect into N partitions
        print 'loading input reads from', inp_name
        ht.consume_fasta(inp_name)

        # Change 0.2 only if you really grok it.  HINT: You don't.
        fp_rate = khmer.calc_expected_collisions(ht)
        print 'fp rate estimated to be %1.3f' % fp_rate

        if fp_rate > 0.20:
            print >>sys.stderr, "**"
            print >>sys.stderr, "** ERROR: the counting hash is too small for"
            print >>sys.stderr, "** %s.  Increase hashsize/num ht." %(inp_name)
            print >>sys.stderr, "**"
            print >>sys.stderr, "** Do not use these results!!"
            print >>sys.stderr, "%s is not processed, inscrease mem" %(inp_name)
            print >> logfile, "%s is not processed, inscrease mem" %(inp_name)
            cnt += 1
            outfp.close()
            os.remove('%s.sweep3' %os.path.basename(inp_name))

        else:
            new_lis.append(query_list[n])

    print '%d files do not have enough mem assigned' %cnt

    print 'starting sweep.'

    n = 0
    m = 0
    for n, record in enumerate(screed.open(readsfile)):
        if len(record.sequence) < K:
            continue

        if n % 100000 == 0:
            print '...', n, m

        for ht, outfp in new_lis:
            count = ht.get_median_count(record.sequence)[0]
            if count:
                outfp.write('>%s\n%s\n' % (record.name, record.sequence))

if __name__ == '__main__':
    main()
