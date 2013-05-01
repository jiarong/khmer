#! /usr/bin/env bash
#PBS -l nodes=1:ppn=1,walltime=168:00:00,mem=1000gb
#PBS -M guojiaro@gmail.com
#PBS -j oe
#PBS -m abe
#PBS -A ged-intel11

###==============> change MEM
# change relative #PBS
MEM_GB=1000
RATIO=0.1 # part of MEM for bashTables
THREADS=1
###==============> repalce the file name

LUMPFILE=/mnt/lustre_scratch_2012/tg/g/data/amo/allIn/P.R1/amo.part.group0515.fa
OUTDIR=$LUMPFILE.deLump.out
SYMLINK=lump.fa
echo

LUMPLABEL=lump
PARTLABEL=amoF.lumpFilt          #PARTLABEL

###===============> REQUEST THIS memory
# keep 4/5 of the MEM for hashTables
MEM=$(bc <<< "$MEM_GB*10^9") ### gb
PART_HASHSIZE=$(bc <<< "scale=2; $MEM*$RATIO/4*8")

echo "$MEM"
echo "$PART_HASHSIZE"
echo "$LUMPLABEL"
echo "$PARTLABEL"

set -e
module load screed

if [ ! -d $OUTDIR ]; then
  mkdir $OUTDIR
else
  echo "$OUTDIR already exists.."
fi
cd $OUTDIR

if [ -h $SYMLINK ]; then
  rm $SYMLINK
fi

ln -s $LUMPFILE $SYMLINK

#Initial round
echo "Initial round (load-graph.py):"
# create graph,
time(
python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/load-graph.py -k 32 -N 4 -x $PART_HASHSIZE $LUMPLABEL $SYMLINK
)
# create an initial set of stoptags to help in knot-traversal; otherwise,
# partitioning and knot-traversal (which is systematic) is really expensive.
python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/make-initial-stoptags.py $LUMPLABEL

echo "Partition graph (partition-graph.py) using the stoptags file"
# now partition the graph, using the stoptags file
time(
python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/partition-graph.py --threads $THREADS --stoptags $LUMPLABEL.stoptags $LUMPLABEL
)

echo "use the partitioned subsets to find the k-mers that nucleate the lump"
# use the partitioned subsets to find the k-mers that nucleate the lump
python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/find-knots.py -x $PART_HASHSIZE -N 4 $LUMPLABEL

echo "remove those k-mers from the fasta files"
# remove those k-mers from the fasta files
python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/filter-stoptags.py *.stoptags $SYMLINK

# now, reload the filtered data set in and partition again.
echo "now, reload the filtered data set in and partition again"
python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/load-graph.py -x $PART_HASHSIZE -N 4 $PARTLABEL $SYMLINK.stopfilt
python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/partition-graph.py --threads $THREADS $PARTLABEL
python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/merge-partitions.py $PARTLABEL
python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/annotate-partitions.py $PARTLABEL $SYMLINK.stopfilt
python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/extract-partitions.py $PARTLABEL $SYMLINK.stopfilt.part
