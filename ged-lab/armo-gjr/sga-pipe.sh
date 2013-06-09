#! /usr/bin/env bash

###
# usage: bash <thisFile> <group####.fa> <overlapCutOff> <outDir>
#
set -e

# Set SEQ to be the path of seq file
#

if [ ! -f $1 ]
then
  echo '$SEQ is not found..'
  exit
fi
SEQ=$(readlink -f $1)
ASSEMBLE_OVERLAP=$2
OUTDIR=$3
MIN_OVERLAP=$((ASSEMBLE_OVERLAP-10))


#
# The number of threads to use
#
CPU=1
ERRATE=0.01

fName=$(basename $SEQ)

PREFIX=$fName.sga.$ASSEMBLE_OVERLAP


if [ ! -d $OUTDIR ]
then
    echo "$OUTDIR doesnot exist.."
    mkdir "$OUTDIR"
    echo "mkdir $OUTDIR"
fi
echo "cd $OUTDIR"
cd $OUTDIR

#
# Parameters
#

# Program paths
SGA_BIN=sga
BWA_BIN=bwa
SAMTOOLS_BIN=samtools
BAM2DE_BIN=sga-bam2de.pl
ASTAT_BIN=sga-astat.py
DISTANCE_EST=DistanceEst # command in ABySS package


# The minimum overlap to use when computing the graph.
# The final assembly can be performed with this overlap or greater
#MIN_OVERLAP=85
#MIN_OVERLAP=$2

# The overlap value to use for the final assembly
#ASSEMBLE_OVERLAP=111
#ASSEMBLE_OVERLAP=$((MIN_OVERLAP+10))

# Branch trim length
TRIM_LENGTH=400

# The minimum length of contigs to include in a scaffold
MIN_CONTIG_LENGTH=400

# The minimum number of reads pairs required to link two contigs
#MIN_PAIRS=10
MIN_PAIRS=4

#
# load tools required
#

module load bwa
module load SAMTools
module load screed
#module load ABySS  # need abyss-fixmate, DistanceEst # default module does not work
module load ABySS/1.3.3
source /mnt/home/guojiaro/Documents/vEnv/sga/bin/activate

echo "###This is the log file for this run, mainly to keep track of running time" > $PREFIX.runLog
#
# Dependency checks
#

# Check the required programs are installed and executable
prog_list="$SGA_BIN $BWA_BIN $SAMTOOLS_BIN $BAM2DE_BIN $DISTANCE_EST $ASTAT_BIN"
for prog in $prog_list; do
    hash $prog 2>/dev/null || { echo "Error $prog not found. Please place $prog on your PATH or update the *_BIN variables in this script"; exit 1; }
done 

# Check the files are found
#file_list="$IN1 $IN2"
file_list="$SEQ"
for input in $file_list; do
    if [ ! -f $input ]; then
        echo "Error input file: $input not found"; exit 1;
    fi
done

# No preprocessing and error correction needed in pipe

#
# Primary (contig) assembly
#

if [ -f $PREFIX.runLog ]; then
  rm $PREFIX.runLog
fi

# Index the corrected data.
echo "$SGA_BIN index -a ropebwt -t $CPU -p $PREFIX $SEQ" >> $PREFIX.runLog
(/usr/bin/time -f "eTime: %E; CPU: %P; maxMem: %M; AveMem: %K" \
$SGA_BIN index -a ropebwt -t $CPU -p $PREFIX $SEQ
)

# Remove exact-match duplicates and reads with low-frequency k-mers

#
#ATT: turn off kmer abundance filter by -k 1 -x 1
#
echo "$SGA_BIN filter -k 1 -x 1 -t $CPU -p $PREFIX --homopolymer-check --low-complexity-check $SEQ" >> $PREFIX.runLog
(/usr/bin/time -f "eTime: %E; CPU: %P; maxMem: %M; AveMem: %K" \
$SGA_BIN filter -k 1 -x 1 -t $CPU -p $PREFIX --homopolymer-check --low-complexity-check $SEQ
)

# Compute the structure of the string graph
echo "$SGA_BIN overlap -m $MIN_OVERLAP -e $ERRATE -t $CPU $PREFIX.filter.pass.fa" >> $PREFIX.runLog
(/usr/bin/time -f "eTime: %E; cpu: %P; maxMem: %M; aveMem: %K" \
$SGA_BIN overlap -m $MIN_OVERLAP -e $ERRATE -t $CPU $PREFIX.filter.pass.fa
)

# Perform the contig assembly
echo "$SGA_BIN assemble -m $ASSEMBLE_OVERLAP --min-branch-length $TRIM_LENGTH -o $PREFIX $PREFIX.filter.pass.asqg.gz" >> $PREFIX.runLog
(/usr/bin/time -f "eTime: %E; cpu: %P; maxMem: %M; aveMem: %K" \
$SGA_BIN assemble -m $ASSEMBLE_OVERLAP --min-branch-length $TRIM_LENGTH -d 0.1 -g 0.05 -o $PREFIX $PREFIX.filter.pass.asqg.gz
)

### job completed

# output job stats
if [ -n "${PBS_JOBID}" ]
then
    qstat -f "${PBS_JOBID}"
fi

echo "job completed.."
