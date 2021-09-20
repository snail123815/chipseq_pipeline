#!/bin/bash
#SBATCH --job-name=bamcoverage
#SBATCH --output=slurm.coverage.out
#SBATCH --error=slurm.coverage.err
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --partition=cpu-medium
#SBATCH --time=06:00:00
#SBATCH --mail-user="c.du@biology.leidenuniv.nl"
#SBATCH --mail-type="ALL"

SOURCEDIR=/home/duc/data/sco1447_chipseq
INDIR=alignmentBAM
OUTFILE=coverage.tsv
SCRIPTDIR=/data/projects/pi-vriesendorpb/duc/chipseq_pipeline

# Make conda avaliable. If you don't have conda installed in your
# own dir, try use "module load Miniconda3"
BASHRC=/home/duc/.bashrc


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
conda activate /data/projects/pi-vriesendorpb/tools/anaconda3/envs/shortReads

# SCRATCH = fast IO location in /scratchdata
# Copy between Scratch-shared to Local scratch (SCRATCH) is very fast.
# (seconds for 1 GB)
SCRATCH=/scratchdata/${SLURM_JOB_USER}/${SLURM_JOB_ID} # fresh for every job
echo "[$SHELL] Node scratch: "$SCRATCH
export RUNDIR=$SCRATCH/chipanalysis
mkdir $RUNDIR
echo "[$SHELL] Run directory: "$RUNDIR

echo "[$SHELL] Copy file to run directory."
cp -r $SOURCEDIR/$INDIR $RUNDIR/
cp $SCRIPTDIR/calculate_coverage.py $RUNDIR/
cd $RUNDIR
echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - $START_TIME))" +%T) (HH:MM:SS)"

echo "[$SHELL] calculating coverage"
echo "[$SHELL] python calculate_coverage.py -p $INDIR -o $OUTFILE" 
python calculate_coverage.py -p $INDIR -o $OUTFILE

echo "[$SHELL] done"
echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - $START_TIME))" +%T) (HH:MM:SS)"
echo "[$SHELL] copy result back"

cp -r $RUNDIR/$OUTFILE $SOURCEDIR/
echo "[$SHELL] done"
echo "[$SHELL] Time elaspsed: $(date -ud "@$(($(date +%s) - $START_TIME))" +%T) (HH:MM:SS)"

