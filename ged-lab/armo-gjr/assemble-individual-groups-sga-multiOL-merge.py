#! /usr/bin/env python
import sys
import tempfile
import screed
from screed.fasta import fasta_iter
import shutil
import os.path
import subprocess
import glob

'''
usage: %%python <thisFile> <group.fasta>
'''

LENGTH_CUTOFF=400
SGA_PIPE='/mnt/home/guojiaro/Documents/lib/git/khmer/ged-lab/armo-gjr/sga-pipe.sh'
print >> sys.stderr, '### make sure sga-pipe.sh path is right: %s' %SGA_PIPE
scripts_dir = os.path.dirname(__file__)
scripts_dir = os.path.abspath(scripts_dir)

def sep_PEandSE(f):
    # produce .pe and .se files
    p = subprocess.Popen('python %s/strip-and-split-for-assembly-gjr.py %s' %(scripts_dir, f), shell=True)
    p.communicate()
    assert p.returncode == 0

def assemble_sequences(f, k, dir, length_cutoff=LENGTH_CUTOFF):
    TAG = os.path.basename(f)
    seqfile = f
    assemble_dir = dir

    p = subprocess.Popen('bash %s %s %d %s' % (SGA_PIPE, seqfile, k, assemble_dir), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #print 'bash %s %s %d %s' % (SGA_PIPE, seqfile, k, assemble_dir)
    (stdout, stderr) = p.communicate()
    assert p.returncode == 0, (stdout, stderr)

    total = 0
    contig_file = os.path.join(assemble_dir, '%s.sga.%d-contigs.fa' %(TAG, k))
    for r in screed.open(contig_file):
        seqlen = len(r['sequence'])
        if seqlen >= length_cutoff:
            total += seqlen

    return total

def best_merge_assemble_sequences(f, assemble_dir, try_k=(35, 40, 50, 70,)):

    #sep_PEandSE(f) # produce .pe and se files
    TAG = os.path.basename(f)
    dest_dir = os.path.dirname(f)

    best_k = try_k[0]
    best_total = assemble_sequences(f, best_k, assemble_dir)
    best_asse ='%s/%s.sga.%d-contigs.fa.best' %(assemble_dir, TAG, best_k)
    os.rename('%s/%s.sga.%d-contigs.fa' %(assemble_dir,TAG,best_k), best_asse)
    print 'total: %.2f Mbp\tk: %d\t%s.sga.%d-contigs.fa' %(float(best_total)/1e6, best_k, TAG, best_k)
    
    for k in try_k[1:]:
        total = assemble_sequences(f, k, assemble_dir)
        print 'total: %.2f Mbp\tk: %d\t%s.sga.%d-contigs.fa' %(float(best_total)/1e6, k, TAG, k)

        if total > best_total:
            # <<
            #os.remove('%s.%d.contig.fa' %(f,best_k))
            # << if uncomment os.remove, comment next 2 lines
            old_best_asse ='%s/%s.sga.%d-contigs.fa.best' %(assemble_dir, TAG, best_k)
            os.rename(old_best_asse,'%s/%s.sga.%d-contigs.fa' %(assemble_dir, TAG, best_k))


            best_total = total
            best_k = k
            
            # << if uncomment os.remove, comment next 2 lines
            best_asse ='%s/%s.sga.%d-contigs.fa.best' %(assemble_dir, TAG, k)
            os.rename('%s/%s.sga.%d-contigs.fa' %(assemble_dir, TAG, k), best_asse)

    # merge contigs from diff min Overlap length

    fLis = glob.glob('%s/%s.sga.*-contigs.fa*' %(assemble_dir, TAG))
    assert len(try_k) == len(fLis), 'number of k should be the same as number'\
                                    ' of contig files'
    combined_asse = '%s/%s.sga.multiKCat' %(assemble_dir, TAG)
    with open(combined_asse, 'wb') as fw:
        cnt = 0
        for ff in fLis:
            for rec in screed.open(ff):
                seq = rec.sequence
                print >> fw, '>%d\n%s' %(cnt, seq)
                cnt += 1

    # sga pipe with k = 100 
    merge_k = 100
    merged_total = assemble_sequences(combined_asse, merge_k, assemble_dir)
    merged_asse = '%s/%s.sga.multiKCat.sga.%d-contigs.fa' %(assemble_dir, TAG, merge_k)

    # move cotigs to same dir as group file
    os.rename(merged_asse, '%s/%s.sga-contigs.fa.merged' %(dest_dir, TAG))
    for ff in fLis:
        os.rename(ff, '%s/%s' %(dest_dir, os.path.basename(ff)))

    return best_k, best_total, merged_total

def main():
    f = os.path.abspath(sys.argv[1])
    dest_dir = os.path.dirname(f)

    try:
        assemble_dir = tempfile.mkdtemp(dir=dest_dir)
        assemble_dir = os.path.abspath(assemble_dir)

        group = 'nogroup'
        fName = os.path.basename(f)
        for i in fName.split('.'):
            if 'group' in i:
                group = i

        k, total, merged_total = best_merge_assemble_sequences(f, assemble_dir)
        print
        print 'best assembly for %s: k=%d, %.2f Mbp' % (fName, k, total*1.0/1e6)
        print 'merged assembly for %s: %.2f Mbp' %(fName, merged_total*1.0/1e6)

    finally:
        shutil.rmtree(assemble_dir)


if __name__ == '__main__':
    main()
