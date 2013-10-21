#!/usr/bin/env python
# modified based on oasis_pipeline.py


import sys
import subprocess
import optparse
import argparse
import shutil
import os, errno
import tempfile

sga_bin = 'sga'
cdhit_bin = 'cd-hit-est'
toAmos_bin = 'toAmos'
minimus2_bin = 'minimus2'

##########################################
## Options and defaults
##########################################
def getOptions():
    #parser = optparse.OptionParser('usage: %prog [options] --data "velveth file descriptors"')
    parser = argparse.ArgumentParser(usage='usage: %prog [options] --data "file.fasta"')
    parser.add_argument('-d', '--data',dest='data',help='SGA file descriptors',metavar='FILE', default='')
    parser.add_argument('-o', '--outdir',dest='outdir',help='output directory for temp files',metavar='DIR', default='')
    parser.add_argument('-l', '--minlen',dest='minlen',type=int,help='min length cutoff',default=300)
    parser.add_argument('-r', '--method',dest='method',help='method for merging assemblies',default='sga', choices = ['amos', 'sga'])
    parser.add_argument('-T', '--thread',dest='thread',type=int,help='number of threads to use for cdhit or sga',default=1)
    parser.add_argument('-O', '--overlap',dest='overlap',type=int,help='overlap for merge by sga fm-merge or minimus2',default=40)
    options = parser.parse_args()
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

def renameAndCat(options):
    # combine and rename due to possible same names
    combined_file = '%s/contigs.cat.fa' %(options.outdir)
    with open(combined_file, 'wb') as fw:
        with open(options.data, 'rb') as fp:
            cnt = 0 # name as index
            seq = ''
            for line in fp:
                line = line.rstrip()
                if line.startswith('>'):
                    if len(seq) < options.minlen:
                        seq = ''
                        continue
                    if 'N' in seq:
                        seq = seq.replace('N', 'G')
                        print >> sys.stderr, '%d has Ns in it, replaced'\
                                                         'with Gs' %cnt
                    print >> fw, '>%d\n%s' %(cnt, seq)
                    cnt += 1
                    seq = ''
                else:
                    seq = '%s%s' %(seq,line)

            if len(seq) >= options.minlen:
                print >> fw, '>%d\n%s' %(cnt, seq)
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

def mergeAsseWithAMOS(options):
    # require cd-hit and AMOS
    p = subprocess.Popen([cdhit_bin, '-i', '%s/contigs.cat.fa' %options.outdir, '-o', '%s/contigs.cat.fa.cdhit' %options.outdir, '-c', '1', '-M', '0', '-T', '%d' %options.thread ], stdout=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "cd-hit-est failed\n"

    p = subprocess.Popen([toAmos_bin, '-s', '%s/contigs.cat.fa.cdhit' %(options.outdir), '-o', '%s/contigs.cat.fa.cdhit.afg' %options.outdir], stdout=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "toAmos failed\n"

    p = subprocess.Popen([minimus2_bin, '%s/contigs.cat.fa.cdhit' %options.outdir, '-D', 'MINID=100', '-D', 'OVERLAP=%d' %options.overlap], stdout=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "minumus2 failed\n"

    files = ['%s/contigs.cat.fa.cdhit.fasta' %(options.outdir), '%s/contigs.cat.fa.cdhit.singletons.seq' %(options.outdir)]
    merged_file = '%s/contigs.fa.amos.merged' %options.outdir
    with open(merged_file, 'wb') as merged_fp:
        for f in files:
            shutil.copyfileobj(open(f, 'rb'), merged_fp)

    # also save the combined and cdhit dereplicated seqs
    cat_cdhit_contigs = '%s.cdhit' %(options.data)
    shutil.copy('%s/contigs.cat.fa.cdhit' %(options.outdir), cat_cdhit_contigs)

    # count total bp
    merged_contigs = '%s.amos.merged' %(f1)
    shutil.move(merged_file, merged_contigs)
    merged_total = countTotal(merged_contigs, options.minlen)

    print 'amos merged\ttotal: %.3f' %(float(merged_total)/1e6)
    print 'amos merged assembly for %s: %d bp' % (os.path.basename(options.data), merged_total)


def mergeAsseWithSGA(options):
    cur_dir = os.path.abspath(os.curdir)
    os.chdir(options.outdir)
    #
    # make sure to chdir back to curdir
    #
    p = subprocess.Popen([sga_bin, 'index', '-a', 'sais', '-d', '500000', '-t', '%d' %options.thread, 'contigs.cat.fa'], stdout=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "sga index failed\n"

    p = subprocess.Popen([sga_bin, 'filter', '-t', '%d' %options.thread, '--no-kmer-check', 'contigs.cat.fa'], stdout=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "sga filter failed\n"

    p = subprocess.Popen([sga_bin, 'fm-merge', '-m', '%d' %options.overlap, '-t', '%d' %options.thread, 'contigs.cat.filter.pass.fa'], stdout=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "sga fm-merge failed\n"

    p = subprocess.Popen([sga_bin, 'index', '-a', 'sais', '-d', '500000', '-t', '%d' %options.thread, 'contigs.cat.filter.pass.merged.fa'], stdout=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "sga index after fm-merge failed\n"

    p = subprocess.Popen([sga_bin, 'filter', '-t', '%d' %options.thread, '--no-kmer-check', 'contigs.cat.filter.pass.merged.fa'], stdout=subprocess.PIPE)
    output = p.communicate()
    assert p.returncode == 0, output[0] + "sga filter after fm -merge failed\n"

    #
    os.chdir(cur_dir)
    #

    # count total bp
    merged_total = 0
    merged_contigs = '%s.sga.merged' %(options.data)
    temp_merged_contigs = '%s/contigs.cat.filter.pass.merged.filter.pass.fa' %(options.outdir)
    shutil.move(temp_merged_contigs, merged_contigs)
    merged_total = countTotal(merged_contigs, options.minlen)
    print 'sga fm-merged\ttotal: %.3f' %(float(merged_total)/1e6)
    print 'sga fm-merged assembly for %s: %d bp' % (os.path.basename(options.data), merged_total)

##########################################
## Checking dependencies
##########################################
def checkSGA():
    try:
        p = subprocess.Popen([sga_bin], stdout=subprocess.PIPE)
    except OSError:
        print "Could not find Velvet"
        print "Make sure that it is properly installed on your path"
        sys.exit(1)
    for line in p.stdout:
        items = line.strip().split(' ')
        if items[0] == 'Version:':
            items2 = map(int, items[1].split('.'))
            assert items2 >= [0,9,35], "Velvet must have version 0.9.35 or higher (currently %s)" % items[1]
            return
    assert False


def checkCDHIT():
    try:
        p = subprocess.Popen([cdhit_bin], stdout=subprocess.PIPE)
    except OSError:
        print "Could not find cd-hit-est"
        print "Make sure that it is properly installed on your path"
        sys.exit(1)

def checkAMOS():
    try:
        p = subprocess.Popen([toAmos_bin], stdout=subprocess.PIPE)
    except OSError:
        print "Could not find toAmos"
        print "Make sure that it is properly installed on your path"
        sys.exit(1)

    try:
        p = subprocess.Popen([minimus2_bin], stdout=subprocess.PIPE)
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

    f1 = os.path.abspath(options.data)
    dest_dir = os.path.dirname(f1)
    # make seq file abspath
    options.data = f1
    if not os.path.isfile(f1):
        print '%s does not exist, exiting..'
        sys.exit(1)
    # use first input seq file's dir as dest dir where contigs are 
    outdir = True
    try:
        if not options.outdir:
            outdir = False
            options.outdir = tempfile.mkdtemp(suffix='.sgaout', dir=dest_dir)
        else:
            try:
                os.makedirs(options.outdir)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(options.outdir):
                    pass
                else:
                    raise

        checkSGA()

        renameAndCat(options)

        if options.method == 'amos':
            checkCDHIT
            checkAMOS
            mergeAsseWithAMOS(options)
        elif options.method == 'sga':
            checkSGA
            mergeAsseWithSGA(options)

    finally:
        if not outdir:
            print 'deleting temp_dir: %s' %options.outdir
            shutil.rmtree(options.outdir)
        else:
            print 'intermediate files kept in: %s' %options.outdir
            pass

if __name__ == "__main__":
    main()
