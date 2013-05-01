#! /usr/bin/env python
import sys
import tempfile
from screed.fasta import fasta_iter
import shutil
import os.path
import subprocess

LENGTH_CUTOFF=400
scripts_dir = os.path.dirname(__file__)
scripts_dir = os.path.abspath(scripts_dir)


def sep_PEandSE(f):
    # produce .pe and .se files
    p = subprocess.Popen('python %s/strip-and-split-for-assembly-gjr.py %s' %(scripts_dir, f), shell=True)
    p.communicate()
    assert p.returncode == 0

def assemble_sequences(f, k, length_cutoff=LENGTH_CUTOFF):
    try:
        seqfile = f
        #dirname = os.path.dirname(os.path.abspath(f))
        dirname = tempfile.mkdtemp()

        assemble_dir = os.path.join(dirname, 'assemble')
        p = subprocess.Popen('velveth %s %d -shortPaired %s.pe -short %s.se' % (assemble_dir, k, seqfile, seqfile), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        assert p.returncode == 0, (stdout, stderr)

        p = subprocess.Popen('velvetg %s -read_trkg yes -exp_cov auto -cov_cutoff 0' % (assemble_dir,),
                             shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        assert p.returncode == 0, (stdout, stderr)

        x = []
        total = 0
        for r in fasta_iter(open(os.path.join(assemble_dir, 'contigs.fa'))):
            seqlen = len(r['sequence'])
            if seqlen >= length_cutoff:
                x.append(r)
                total += seqlen

        return total, x
    finally:
        pass
        shutil.rmtree(dirname)
        #print 'XXX', dirname

def best_assemble_sequences(f, try_k=(33, 37, 39, 49, 69)):

    sep_PEandSE(f) # produce .pe and se files

    best_k = try_k[0]
    best_total, best_records = assemble_sequences(f, best_k)
    print 'total: %.2f\tk: %d\t%s' %(float(best_total)/1e6, best_k, os.path.basename(f))
    
    for k in try_k[1:]:
        total, records = assemble_sequences(f, k)
        print 'total: %.2f\tk: %d\t%s' %(float(total)/1e6, k, os.path.basename(f))

        if total > best_total:
            best_total = total
            best_records = records
            best_k = k

    return best_k, best_total, best_records

group = None
for i in  sys.argv[1].split('.'):
    if 'group' in i:
        group = i

k, total, records = best_assemble_sequences(sys.argv[1])
print
print 'best assembly for %s: k=%d, %d bp' % (sys.argv[1], k, total)

fp = open(sys.argv[1] + '.%d.best' %k, 'wb')
for n,r in enumerate(records):
    fp.write('>%s.%d\n%s\n' % (group, n, r['sequence']))

fp.close()
