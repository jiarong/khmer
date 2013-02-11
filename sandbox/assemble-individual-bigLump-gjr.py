#! /usr/bin/env python
import sys
import tempfile
import screed
from screed.fasta import fasta_iter
import shutil
import os.path
import subprocess

LENGTH_CUTOFF=400
scripts_dir = os.path.dirname(__file__)
scripts_dir = os.path.abspath(scripts_dir)

f = os.path.abspath(sys.argv[1])
dest_dir = os.path.dirname(f)


def sep_PEandSE(f):
    # produce .pe and .se files
    p = subprocess.Popen('python %s/strip-and-split-for-assembly-gjr.py %s' %(scripts_dir, f), shell=True)
    p.communicate()
    assert p.returncode == 0

def assemble_sequences(f, k, length_cutoff=LENGTH_CUTOFF):
    try:
        seqfile = f
        #dirname = os.path.dirname(os.path.abspath(f))
        dirname = tempfile.mkdtemp(dir=dest_dir)

        assemble_dir = os.path.join(dirname, 'assemble')
        p = subprocess.Popen('velveth %s %d -shortPaired %s.pe -short %s.se' % (assemble_dir, k, seqfile, seqfile), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        assert p.returncode == 0, (stdout, stderr)

        p = subprocess.Popen('velvetg %s -read_trkg yes -exp_cov auto -cov_cutoff 0' % (assemble_dir,),
                             shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = p.communicate()
        assert p.returncode == 0, (stdout, stderr)

        total = 0
        contig_file = os.path.join(assemble_dir, 'contigs.fa')
        for r in screed.open(contig_file):
            seqlen = len(r['sequence'])
            if seqlen >= length_cutoff:
                total += seqlen

        shutil.move(contig_file, '%s.%d.contig.fa' %(f,k))

        return total

    finally:
        shutil.rmtree(dirname)
        #print 'XXX', dirname

#def best_assemble_sequences(f, try_k=(33, 37, 39, 49, 69)):
def best_assemble_sequences(f, try_k=(49, 69)):

    sep_PEandSE(f) # produce .pe and se files

    best_k = try_k[0]
    best_total = assemble_sequences(f, best_k)
    best_asse ='%s.%d.contig.fa.best' %(f,best_k)
    os.rename('%s.%d.contig.fa' %(f,best_k), best_asse)
    print 'total: %.2f\tk: %d\t%s.%d.contig.fa' %(float(best_total)/1e6, best_k, os.path.basename(f), best_k)
    
    for k in try_k[1:]:
        total = assemble_sequences(f, k)
        print 'total: %.2f\tk: %d\t%s.%d.contig.fa' %(float(best_total)/1e6, k, os.path.basename(f), k)

        if total > best_total:
            # <<
            #os.remove('%s.%d.contig.fa' %(f,best_k))
            # << if uncomment os.remove, comment next 2 lines
            old_best_asse ='%s.%d.contig.fa.best' %(f,best_k)
            os.rename(old_best_asse,'%s.%d.contig.fa' %(f,best_k))


            best_total = total
            best_k = k
            
            # << if uncomment os.remove, comment next 2 lines
            best_asse ='%s.%d.contig.fa.best' %(f,k)
            os.rename('%s.%d.contig.fa' %(f,k), best_asse)


    return best_k, best_total

group = None
fname = os.path.basename(sys.argv[1])
for i in fname.split('.'):
    if 'group' in i:
        group = i

k, total = best_assemble_sequences(f)
print
print 'best assembly for %s: k=%d, %d bp' % (sys.argv[1], k, total)
