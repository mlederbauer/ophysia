#!/bin/bash
#SBATCH --job-name=CCO
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --output=CCO_%j.out
#SBATCH --error=CCO_%j.err

# Create temporary directory in scratch
SCRATCH_WORK="$SCRATCH_BASE/CCO/$SLURM_JOB_ID"
mkdir -p "$SCRATCH_WORK"

# Copy input files to scratch
cp /cluster/home/mlederbauer/cco_analysis/CCO/molecule.inp "$SCRATCH_WORK/"
cp /cluster/home/mlederbauer/cco_analysis/CCO/molecule.xyz "$SCRATCH_WORK/"

# Change to scratch directory
cd "$SCRATCH_WORK"

# Run ORCA
/cluster/software/commercial/orca/6.0.0/x86_64/bin/orca molecule.inp > molecule.out

# Copy results back to home directory
cp *.gbw *.xyz *.out *.hess /cluster/home/mlederbauer/cco_analysis/CCO/

# Cleanup scratch
rm -rf "$SCRATCH_WORK"
