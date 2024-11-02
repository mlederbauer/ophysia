#!/bin/bash
#SBATCH --job-name=CC_cation
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --output=CC_cation_%j.out
#SBATCH --error=CC_cation_%j.err

# Create temporary directory in scratch
SCRATCH_WORK="$SCRATCH_BASE/CC_cation/$SLURM_JOB_ID"
mkdir -p "$SCRATCH_WORK"

# Copy input files to scratch
cp /cluster/home/mlederbauer/cco_analysis/CC_cation/molecule.inp "$SCRATCH_WORK/"
cp /cluster/home/mlederbauer/cco_analysis/CC_cation/molecule.xyz "$SCRATCH_WORK/"

# Change to scratch directory
cd "$SCRATCH_WORK"

# Run ORCA
/cluster/software/commercial/orca/6.0.0/x86_64/bin/orca molecule.inp > molecule.out

# Copy results back to home directory
cp *.gbw *.xyz *.out /cluster/home/mlederbauer/cco_analysis/CC_cation/

# Cleanup scratch
rm -rf "$SCRATCH_WORK"
