#!/bin/bash
#SBATCH --job-name=align
#SBATCH --output=align.out
#SBATCH --error=align.err
#SBATCH --ntasks=16
#SBATCH --mem=8G
#SBATCH --partition=cpu-short
#SBATCH --time=03:00:00

# ============================================================================================
CONDA=micromamba # use CONDA=conda when you have your conda directly in your .bashrc.
 		 # directly using "alias conda=micromamba" in .bashrc do not work with testing Node
		 # maybe because the testing Node is running as combined node by default.
		 # eg: nodelogin[01-02]
		 # setting #SBATCH --nodes=1 should address this problem
NCPU=16
SOURCEDIR=/path/to/sourcedir/with/rawreads/dir/and/genome/file
RAWREADS=raw_reads # dir in source dir
GENOME=M145.fa # file in source dir
ISPE="-pe"

SCRIPTDIR=/data/projects/pi-vriesendorpb/duc/chipseq_pipeline

# Make conda avaliable. If you don't have conda installed in your
# own dir, try use "module load Miniconda3"
BASHRC=/home/duc/.bashrc

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


# Enter the designed environment
echo "[$SHELL] Activate shortReads conda env"
source $BASHRC
$CONDA activate /data/projects/pi-vriesendorpb/tools/anaconda3/envs/shortReads

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

