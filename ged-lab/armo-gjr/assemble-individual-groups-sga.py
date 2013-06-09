#! /usr/bin/env python
import sys
import tempfile
from screed.fasta import fasta_iter
import shutil
import os.path
import subprocess

LENGTH_CUTOFF=400

SGA_PIPE='/mnt/home/guojiaro/Documents/lib/git/khmer/ged-lab/armo-gjr/sga-pipe.sh'
print >> sys.stderr, '### make sure sga-pipe.sh path is right: $SGA_PIPE'

scripts_dir = os.path.dirname(__file__)
scripts_dir = os.path.abspath(scripts_dir)

def assemble_sequences(f, k, length_cutoff=LENGTH_CUTOFF):
    try:
        seqfile = f
        #dirname = os.path.dirname(os.path.abspath(f))
        dirname = tempfile.mkdtemp()

        assemble_dir = os.path.join(dirname, 'assemble')
        p = subprocess.Popen('bash %s %s %d %s' % (SGA_PIPE, seqfile, k, assemble_dir), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print 'bash %s %s %d %s' % (SGA_PIPE, seqfile, k, assemble_dir)
        (stdout, stderr) = p.communicate()
        assert p.returncode == 0, (stdout, stderr)

        x = []
        total = 0
        print os.listdir(assemble_dir)
        for r in fasta_iter(open(os.path.join(assemble_dir, '%s.sga.%d-contigs.fa' %(os.path.basename(f), k)))):
            seqlen = len(r['sequence'])
            if seqlen >= length_cutoff:
                x.append(r)
                total += seqlen

        return total, x
    finally:
        pass
        shutil.rmtree(dirname)
        #print 'XXX', dirname

def best_assemble_sequences(f, try_k=(30, 40, 50, 70)):

    best_k = try_k[0]
    best_total, best_records = assemble_sequences(f, best_k)
    print 'total: %.2f(Mbp)\tk: %d\t%s' %(float(best_total)/1e6, best_k, os.path.basename(f))
    
    for k in try_k[1:]:
        total, records = assemble_sequences(f, k)
        print 'total: %.2f(Mbp)\tk: %d\t%s' %(float(total)/1e6, k, os.path.basename(f))

        if total > best_total:
            best_total = total
            best_records = records
            best_k = k

    return best_k, best_total, best_records

group = 'nogroup' 
for i in  sys.argv[1].split('.'):
    if 'group' in i:
        group = i

k, total, records = best_assemble_sequences(sys.argv[1])
print
print 'best assembly for %s: k=%d, %d bp' % (sys.argv[1], k, total)

fp = open(sys.argv[1] + '.sga.%d.best' %k, 'wb')
for n,r in enumerate(records):
    fp.write('>%s.%d\n%s\n' % (group, n, r['sequence']))

fp.close()
