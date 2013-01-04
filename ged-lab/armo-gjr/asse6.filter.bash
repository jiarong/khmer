#! /usr/bin/env bash
#PBS -l nodes=1:ppn=8,walltime=72:00:00,mem=1000gb
#PBS -M guojiaro@gmail.com
#PBS -j oe
#PBS -m abe
#PBS -A ged-intel11

###==============> change MEM
# change relative #PBS
MEM_GB=1000
RATIO=0.998  # part of MEM for bashTables
THREADS=8
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


# PASS3: filter high abund > 50
echo "PASS2: filter-below-abund.py"
echo "PASS2: filter-below-abund.py" 1>&2
time(
python /mnt/home/guojiaro/Documents/lib/git/khmer/sandbox/filter-below-abund.py $HASHTABLE *.keep
)
