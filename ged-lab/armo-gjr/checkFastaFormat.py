#! /usr/bin/env python
# check bad characters in fasta file, espcially miss end of line
# by gjr; 10/19/13

import sys
import screed

def main():
    '''
    python <thisFile> <file.fa>
    '''
    if len(sys.argv) != 2:
        print >> sys.stderr, 'usage: python %s <file.fa>' \
                                %(os.path.basename(sys.argv[0]))
    fastaF = sys.argv[1]

    for n, record in enumerate(screed.open(fastaF)):
        seq = record.sequence.upper()
        if not seq.isalpha():
            print >> sys.stderr, 'Bad character detected at line %d' %(2*n+1)
            print >> sys.stderr, '>%s\n%s' %(repr(record.name), repr(seq))

    print >> sys.stderr, 'Done. %d seqs scanned..' %(n+1)

if __name__ == '__main__':
    main()
