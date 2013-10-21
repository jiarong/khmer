#!/usr/bin/env python
# modified based on oasis_pipeline.py


import sys
import subprocess
import optparse
import shutil
import os, errno
import tempfile

temp_dir = None


##########################################
## Options and defaults
##########################################
def getOptions():
    parser = optparse.OptionParser('usage: %prog [options] --data "velveth file descriptors"')
    parser.add_option('-d', '--data',dest='data',help='Velveth file descriptors',metavar='FILES', default='')
    parser.add_option('-o', '--outdir',dest='outdir',help='output directory for temp files',metavar='DIR', default='')
    parser.add_option('-m', '--kmin',dest='kmin',type="int",help='Minimum k',default=29)
    parser.add_option('-M', '--kmax',dest='kmax',type="int",help='Maximum k',default=69)
    parser.add_option('-s', '--kstep',dest='kstep',type="int",help='Steps in k',default=10)
    parser.add_option('-g', '--merge',dest='kmerge',type="int",help='Merge k',default=29)
    parser.add_option('-l', '--minlen',dest='minlen',type="int",help='min length cutoff',default=300)
    parser.add_option('-r', '--method',dest='method',help='method for merging assemblies',default='sga')
    parser.add_option('-T', '--thread',dest='thread',type="int",help='number of threads to use for cdhit or sga',default=1)
    parser.add_option('-O', '--overlap',dest='overlap',type="int",help='overlap for merge by sga fm-merge or minimus2',default=40)
    options, args = parser.parse_args()
    if  len(options.data) == 0:
        parser.print_help()
        print ''
        print 'You forgot to provide some data files!'
        print 'Current options are:'
        print options
        sys.exit(1)
    return options

##########################################
## Assembly procedure
##########################################
def singleKAssemblies(options):
    #global temp_dir
    p = subprocess.Popen(['velveth', '%s/assembly' %(temp_dir), '%i,%i,%i' % (options.kmin, options.kmax, options.kstep)] + options.data.split(), stdout=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "Hash failed\n"
    for k in range(options.kmin, options.kmax, options.kstep):
        p = subprocess.Popen(['velvetg','%s/assembly_%i' % (temp_dir, k)], stdout=subprocess.PIPE)
        output = p.communicate()
        assert p.returncode == 0, "Velvetg failed at k = %i\n%s" % (k, output[0])

    getBestK(options)


def getBestK(options):
    best_k = 0
    best_total = 0
    for k in range(options.kmin, options.kmax, options.kstep):
        total = 0
        seq = ''
        for line in open('%s/assembly_%d/contigs.fa' %(temp_dir, k)):
            line = line.rstrip()
            if line.startswith('>'):
                if len(seq) < options.minlen:
                    seq = ''
                    continue
                total += len(seq)
                seq = ''
            else:
                seq = '%s%s' %(seq,line)

        if len(seq) >= options.minlen:
            total += len(seq)
            seq = ''

        if total > best_total:
            best_k = k
            best_total = total
        print 'k: %d\ttotal: %.3f' %(k, float(total)/1e6)

    f1 = os.path.abspath(options.data.split()[0])
    print 'best assembly for %s: k=%d, %d bp' % (os.path.basename(f1), best_k, best_total)
    if best_k == 0: 
        print 'nothing assembled with min length %d' %(options.minlen)
        print 'exiting..'
        sys.exit(1)

    # rename and filter short contigs < options.minlen
    bestk_contigs = '%s.%d.vel.contig.fa.best' %(f1, best_k)
    with open(bestk_contigs, 'wb') as fw:
        with open('%s/assembly_%d/contigs.fa' %(temp_dir, best_k), 'rb') as fp:
            cnt = 0 # name as index
            seq = ''
            for line in fp:
                line = line.rstrip()
                if line.startswith('>'):
                    if len(seq) < options.minlen:
                        seq = ''
                        continue
                    print >> fw, '>%d\n%s' %(cnt, seq)
                    cnt += 1
                    seq = ''
                else:
                    seq = '%s%s' %(seq,line)

            if len(seq) >= options.minlen:
                print >> fw, '>%d\n%s' %(cnt, seq)
                cnt += 1
                seq = ''

def renameAndCat(options):
    # combine and rename due to possible same names
    combined_file = '%s/contigs.cat.fa' %(temp_dir)
    with open(combined_file, 'wb') as fw:
        for X in range(options.kmin, options.kmax, options.kstep):
            with open("%s/assembly_%i/contigs.fa" % (temp_dir, X), 'rb') as fp:
                cnt = 0 # name as index
                seq = ''
                for line in fp:
                    line = line.rstrip()
                    if line.startswith('>'):
                        if len(seq) < options.minlen:
                            seq = ''
                            continue
                        print >> fw, '>%d_%d\n%s' %(X, cnt, seq)
                        cnt += 1
                        seq = ''
                    else:
                        seq = '%s%s' %(seq,line)

                if len(seq) >= options.minlen:
                    print >> fw, '>%d_%d\n%s' %(X, cnt, seq)
                    cnt += 1
                    seq = ''

def countTotal(seq, minlen):
    merged_total = 0
    with open(seq, 'rb') as fp:
        cnt = 0 # name as index
        seq = ''
        for line in fp:
            line = line.rstrip()
            if line.startswith('>'):
                if len(seq) < minlen:
                    seq = ''
                    continue
                merged_total += len(seq)
                cnt += 1
                seq = ''
            else:
                seq = '%s%s' %(seq,line)

        if len(seq) >= minlen:
            merged_total += len(seq)
            cnt += 1
            seq = ''

    return merged_total

def mergeAssemblies(options):
    global temp_dir
    p = subprocess.Popen(['velveth','%s/merged' % temp_dir, str(options.kmerge), '-long', '%s/contigs.cat.fa' %temp_dir], stdout=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "Velveth failed at merge\n"
    p = subprocess.Popen(['velvetg','%s/merged' % temp_dir,'-conserveLong','yes', '-min_contig_lgth', '%d' %options.minlen], stdout=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "Velvetg failed at merge\n"
    # mv contig file to '%s/contigs.fa.merged'
    f1 = os.path.abspath(options.data.split()[0])
    shutil.move('%s/merged/contigs.fa' %(temp_dir), '%s.vel.contig.fa.vel.merged' %f1)

def mergeAsseWithAMOS(options):
    # require cd-hit and AMOS
    global temp_dir
    p = subprocess.Popen(['cd-hit-est', '-i', '%s/contigs.cat.fa' %temp_dir, '-o', '%s/contigs.cat.fa.cdhit' %temp_dir, '-c', '1', '-M', '0', '-T', '%d' %options.thread ], stdout=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "cd-hit-est failed\n"

    p = subprocess.Popen(['toAmos', '-s', '%s/contigs.cat.fa.cdhit' %(temp_dir), '-o', '%s/contigs.cat.fa.cdhit.afg' %temp_dir], stdout=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "toAmos failed\n"

    p = subprocess.Popen(['minimus2', '%s/contigs.cat.fa.cdhit' %temp_dir, '-D', 'MINID=100', '-D', 'OVERLAP=%d' %options.overlap], stdout=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "minumus2 failed\n"

    files = ['%s/contigs.cat.fa.cdhit.fasta' %(temp_dir), '%s/contigs.cat.fa.cdhit.singletons.seq' %(temp_dir)]
    merged_file = '%s/contigs.fa.amos.merged' %temp_dir
    with open(merged_file, 'wb') as merged_fp:
        for f in files:
            shutil.copyfileobj(open(f, 'rb'), merged_fp)

    f1 = os.path.abspath(options.data.split()[0])

    # also save the combined and cdhit dereplicated seqs
    cat_cdhit_contigs = '%s.vel.contig.fa.multiKcat.cdhit' %(f1)
    shutil.copy('%s/contigs.cat.fa.cdhit' %(temp_dir), cat_cdhit_contigs)

    # count total bp
    merged_contigs = '%s.vel.contig.fa.amos.merged' %(f1)
    shutil.move(merged_file, merged_contigs)
    merged_total = countTotal(merged_contigs, options.minlen)

    print 'amos merged\ttotal: %.3f' %(float(merged_total)/1e6)
    print 'amos merged assembly for %s: %d bp' % (os.path.basename(f1), merged_total)


def mergeAsseWithSGA(options):
    # require cd-hit and AMOS
    global temp_dir
    cur_dir = os.path.abspath(os.curdir)
    os.chdir(temp_dir)
    #
    # make sure to chdir back to curdir
    #
    p = subprocess.Popen(['sga', 'index', '-a', 'ropebwt', '-t', '%d' %options.thread, 'contigs.cat.fa'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "sga index failed\n"

    p = subprocess.Popen(['sga', 'filter', '-t', '%d' %options.thread, '--no-kmer-check', '--homopolymer-check', '--low-complexity-check', 'contigs.cat.fa'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "sga filter failed\n"

    p = subprocess.Popen(['sga', 'fm-merge', '-m', '%d' %options.overlap, '-t', '%d' %options.thread, 'contigs.cat.filter.pass.fa'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "sga fm-merge failed\n"

    p = subprocess.Popen(['sga', 'index', '-a', 'ropebwt', '-t', '%d' %options.thread, 'contigs.cat.filter.pass.merged.fa'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "sga index after fm-merge failed\n"

    p = subprocess.Popen(['sga', 'filter', '-t', '%d' %options.thread, '--no-kmer-check', '--homopolymer-check', '--low-complexity-check', 'contigs.cat.filter.pass.merged.fa'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "sga filter after fm -merge failed\n"

    #
    os.chdir(cur_dir)
    #

    f1 = os.path.abspath(options.data.split()[0])

    # count total bp
    merged_total = 0
    merged_contigs = '%s.vel.contig.fa.sga.merged' %(f1)
    temp_merged_contigs = '%s/contigs.cat.filter.pass.merged.filter.pass.fa' %(temp_dir)
    shutil.move(temp_merged_contigs, merged_contigs)
    merged_total = countTotal(merged_contigs, options.minlen)
    print 'sga fm-merged\ttotal: %.3f' %(float(merged_total)/1e6)
    print 'sga fm-merged assembly for %s: %d bp' % (os.path.basename(f1), merged_total)

##########################################
## Checking dependencies
##########################################
def checkVelvet():
    try:
        p = subprocess.Popen(['velveth'], stdout=subprocess.PIPE)
    except OSError:
        print "Could not find Velvet"
        print "Make sure that it is properly installed on your path"
        sys.exit(1)
    for line in p.stdout:
        items = line.strip().split(' ')
        if items[0] == 'Version':
            items2 = map(int, items[1].split('.'))
            assert items2 >= [1,1,7], "Velvet must have version 1.1.07 or higher (currently %s)" % items[1]
            return
    assert False

def checkSGA():
    try:
        p = subprocess.Popen(['sga'], stdout=subprocess.PIPE)
    except OSError:
        print "Could not find SGA"
        print "Make sure that it is properly installed on your path"
        sys.exit(1)

def checkCDHIT():
    try:
        p = subprocess.Popen(['cd-hit-est'], stdout=subprocess.PIPE)
    except OSError:
        print "Could not find cd-hit-est"
        print "Make sure that it is properly installed on your path"
        sys.exit(1)

def checkAMOS():
    try:
        p = subprocess.Popen(['minimus2'], stdout=subprocess.PIPE)
    except OSError:
        print "Could not find minimus2"
        print "Make sure that it is properly installed on your path"
        sys.exit(1)

##########################################
## Clean up
##########################################
def clean(options):
    for k in range(options.kmin, options.kmax, options.kstep):
        shutil.rmtree("%s_%i" % (options.directoryRoot, k))
    os.remove("%sMerged/Sequences" % options.directoryRoot)
    os.remove("%sMerged/Roadmaps" % options.directoryRoot)
    os.remove("%sMerged/PreGraph" % options.directoryRoot)
    os.remove("%sMerged/Graph2" % options.directoryRoot)
    os.remove("%sMerged/LastGraph" % options.directoryRoot)
    os.remove("%sMerged/contigs.fa" % options.directoryRoot)
    os.remove("%sMerged/Log" % options.directoryRoot)

##########################################
## Master function
##########################################
def main():
    options = getOptions()
    assert options.method in ['velvet', 'amos', 'sga'], '%s is not supported' \
                                                         %options.method
    # make kmax inclusive
    if (options.kmax - options.kmin)%options.kstep == 0:
        options.kmax += 1
    
    f1 = os.path.abspath(options.data.split()[0])
    dest_dir = os.path.dirname(f1)
    if not os.path.isfile(f1):
        print '%s does not exist, exiting..'
        sys.exit(1)
    # use first input seq file's dir as dest dir where contigs are 
    try:
        global temp_dir
        if not options.outdir:
            temp_dir = tempfile.mkdtemp(suffix='.velvetout', dir=dest_dir)
        else:
            temp_dir = options.outdir
            try:
                os.makedirs(temp_dir)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(temp_dir):
                    pass
                else:
                    raise

        checkVelvet()
        singleKAssemblies(options)

        renameAndCat(options)
        if options.method == 'velvet':
            mergeAssemblies(options)
        elif options.method == 'amos':
            checkCDHIT
            checkAMOS
            mergeAsseWithAMOS(options)
        elif options.method == 'sga':
            checkSGA
            mergeAsseWithSGA(options)

    finally:
        if not options.outdir:
            print 'deleting temp_dir: %s' %temp_dir
            shutil.rmtree(temp_dir)
        else:
            print 'intermediate files kept in: %s' %temp_dir
            pass

if __name__ == "__main__":
    main()
