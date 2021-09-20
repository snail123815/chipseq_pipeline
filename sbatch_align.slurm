#!/bin/bash
#SBATCH --job-name=align
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --ntasks=16
#SBATCH --mem=8G
#SBATCH --partition=testing
#SBATCH --time=01:00:00
#SBATCH --mail-user="c.du@biology.leidenuniv.nl"
#SBATCH --mail-type="ALL"

NCPU=16
SOURCEDIR=/home/duc/data/sco1447_chipseq
RAWREADS=raw_reads
GENOME=M145.fa
ISPE="-pe"

SCRIPTDIR=/data/projects/pi-vriesendorpb/duc/chipseq_pipeline

# Make conda avaliable. If you don't have conda installed in your
# own dir, try use "module load Miniconda3"
BASHRC=/hoem/duc/.bashrc

OUTDIRSAM=alignmentSAM
OUTDIRBAM=alignmentBAM # this will be actually copied back to SOURCEDIR. 

# ===========================================================================================
# Adding [$SHELL] will give distinction between echo and script output
echo "[$SHELL] #### Starting Job"
# The SLURM_JOB_USER and SLURM_JOB_ID is automatically obtained from sbatch command
echo "[$SHELL] This is $SLURM_JOB_USER, job ID $SLURM_JOB_ID"

# Get start time
START_TIME=$(date +%s)
echo "[$SHELL] Started at: $(date)"

# Change language settings in case perl complains
export LANGUAGE=en_US.UTF-8
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

# Assign CWD, not required for this task
export CWD=$(pwd)
echo "[$SHELL] CWD: "$CWD
source $BASHRC


# Enter the designed environment
echo "[$SHELL] Activate shortReads conda env"
conda activate /data/projects/pi-vriesendorpb/tools/anaconda3/envs/shortReads

# SCRATCH = fast IO location in /scratchdata
# Copy between Scratch-shared to Local scratch (SCRATCH) is very fast.
# (seconds for 1 GB)
SCRATCH=/scratchdata/${SLURM_JOB_USER}/${SLURM_JOB_ID} # fresh for every job
echo "[$SHELL] Node scratch: "$SCRATCH
RUNDIR=$SCRATCH/chipanalysis
mkdir $RUNDIR
echo "[$SHELL] Run directory: "$RUNDIR

echo "[$SHELL] Copy files to run directory."
cp -r $SOURCEDIR/$RAWREADS $RUNDIR/
cp $SOURCEDIR/$GENOME $RUNDIR/
cp $SCRIPTDIR/align_to_genome.py $RUNDIR/
cp $SCRIPTDIR/processSam2Bam.py $RUNDIR/
cd $RUNDIR
echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - $START_TIME))" +%T) (HH:MM:SS)"

echo "[$SHELL] run alignment."
CMD="python align_to_genome.py -r $RAWREADS -g $GENOME -o $OUTDIRSAM -p $NCPU $ISPE"
echo "[$SHELL] $CMD"
eval "$CMD" 

echo "[$SHELL] done"
echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - $START_TIME))" +%T) (HH:MM:SS)"

echo "[$SHELL] run indexing"
CMD="python processSam2Bam.py -p $OUTDIRSAM -o $OUTDIRBAM -t $NCPU"
echo "[$SHELL] $CMD"
eval "$CMD" 

echo "[$SHELL] done"
echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - $START_TIME))" +%T) (HH:MM:SS)"
echo "[$SHELL] copy result back"

cp -r $RUNDIR/$OUTDIRBAM $SOURCEDIR/
echo "[$SHELL] done"
echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - $START_TIME))" +%T) (HH:MM:SS)"

