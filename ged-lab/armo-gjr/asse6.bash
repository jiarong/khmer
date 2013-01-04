#! /usr/bin/env bash
#PBS -l nodes=1:ppn=1,walltime=168:00:00,mem=1000gb
#PBS -M guojiaro@gmail.com
#PBS -j oe
#PBS -m abe
#PBS -A ged-intel11

###==============> change MEM
# change relative #PBS
MEM_GB=1000
RATIO=0.998  # part of MEM for bashTables
THREADS=1
DIGINORM_C=10
###==============> repalce the file name
OUTDIR=/mnt/scratch/tg/g/data/amo/allIn/P.R1
#OUTDIR=/mnt/scratch/tg/g/test
LIS=$(find /mnt/scratch/tg/g/data/amo/P.R1.A* -name *.afterMerge.fastq)
echo "Number of files to process:"
find /mnt/scratch/tg/g/data/amo/P.R1.A* -name *.afterMerge.fastq | wc -l
echo
#LIS=$(find /mnt/scratch/tg/g/test -name 10K*.fastq)

HASHTABLE=hash.hk          #HASHTABLE
PARTLABEL=amo.part          #PARTLABEL

###===============> REQUEST THIS memory
# keep 4/5 of the MEM for hashTables
MEM=$MEM_GB*10^9 ### gb
MEM=$(echo "$MEM" | bc)
DIGNORM_HASHSIZE=$(echo "scale=2; $MEM*$RATIO/4" | bc)
PART_HASHSIZE=$(echo "scale=2; $MEM*$RATIO/4*8" | bc)

echo "$MEM"
echo "$DIGNORM_HASHSIZE"
echo "$PART_HASHSIZE"
echo "$HASHTABLE"
echo "$PARTLABEL"

set -e
module load screed

if [ ! -d $OUTDIR ]; then
  mkdir $OUTDIR
else
  echo "$OUTDIR already exists.."
fi

cd $OUTDIR
CNT=0

# PASS1: digiNorm -C $DIGINORM_C
echo PASS1 digiNorm -C $DIGINORM_C:
echo PASS1 digiNorm -C $DIGINORM_C: 1>&2

for i in $LIS
do
  SEQ=$i

  if [ $CNT -eq 0 ]; then
    time(
    python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/normalize-by-median.py -k 20 -C $DIGINORM_C -x $DIGNORM_HASHSIZE -N 4 -R $HASHTABLE.$CNT.report --savehash $HASHTABLE  $SEQ
    echo $SEQ processed by normalzie-by-median.py > $HASHTABLE.log
    )
  else
    time(
    python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/normalize-by-median.py -k 20 -C $DIGINORM_C -x $DIGNORM_HASHSIZE -N 4 -R $HASHTABLE.$CNT.report --savehash $HASHTABLE -l $HASHTABLE $SEQ
    )
    echo $SEQ processed by normalzie-by-median.py >> $HASHTABLE.log
  fi

  CNT=$((CNT+1))

done

# PASS3: filter high abund > 50
echo "PASS2: filter-below-abund.py"
echo "PASS2: filter-below-abund.py" 1>&2
time(
python /mnt/home/guojiaro/Documents/lib/git/khmer/sandbox/filter-below-abund.py $HASHTABLE *.keep
)
### Partitioning
echo "start partitioning"
echo "start partitioning" 1>&2
#Initial round
echo "Initial round (load-graph.py):"
echo "Initial round (load-graph.py):" 1>&2
time(
python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/load-graph.py -k 32 -N 4 -x $PART_HASHSIZE $PARTLABEL *.keep.below
)
echo "Partition graph (partition-graph.py):"
echo "Partition graph (partition-graph.py):" 1>&2
time(
python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/partition-graph.py --threads 8 -s 1e6 $PARTLABEL
)

echo "Merge parts (merge-partitions.py):"
echo "Merge parts (merge-partitions.py):" 1>&2
time(
python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/merge-partitions.py $PARTLABEL
)

echo "Annotate parts (annotate-partitions.py):"
echo "Annotate parts (annotate-partitions.py):" 1>&2
time(
python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/annotate-partitions.py $PARTLABEL *.keep.below
)

echo "extract-partitions.py:"
echo "extract-partitions.py:" 1>&2
time(
python /mnt/home/guojiaro/Documents/lib/git/khmer/scripts/extract-partitions.py $PARTLABEL *.keep.below.part
)
