#!/bin/bash

# Directory setup
HOME_DIR="/cluster/home/$USER/cco_analysis"
SCRATCH_DIR="/cluster/scratch/$USER/cco_analysis"

# Job parameters
WALLTIME="24:00:00"
MEM_PER_CPU="4000"
NTASKS="24"

# Create necessary directories
mkdir -p $HOME_DIR
mkdir -p $SCRATCH_DIR

# Load required modules
module purge
module load stack/2024-06
module load gcc/12.2.0
module load orca/6.0.0
module load openmpi/4.1.6

# Base directories
HOME_DIR="$PWD/cco_analysis"
SCRATCH_BASE="/cluster/scratch/$USER/cco_analysis"

# Ensure scratch base directory exists
mkdir -p "$SCRATCH_BASE"

# Function to submit a job
submit_orca_job() {
    local dir=$1
    local basename=$(basename "$dir")
    local scratch_dir="$SCRATCH_BASE/$basename"
    
    # Create scratch directory
    mkdir -p "$scratch_dir"
    
    cd "$dir"
    
    if [ ! -f "molecule.inp" ]; then
        echo "No input file found in $dir"
        return
    fi
    
    # Create submission script
    cat > "run_job.sh" << EOF
#!/bin/bash
#SBATCH --job-name=${basename}
#SBATCH --nodes=1
#SBATCH --ntasks=${NTASKS}
#SBATCH --time=${WALLTIME}
#SBATCH --mem-per-cpu=${MEM_PER_CPU}
#SBATCH --output=${basename}_%j.out
#SBATCH --error=${basename}_%j.err

# Create temporary directory in scratch
SCRATCH_WORK="\$SCRATCH_BASE/$basename/\$SLURM_JOB_ID"
mkdir -p "\$SCRATCH_WORK"

# Copy input files to scratch
cp $HOME_DIR/$basename/molecule.inp "\$SCRATCH_WORK/"
cp $HOME_DIR/$basename/molecule.xyz "\$SCRATCH_WORK/"

# Change to scratch directory
cd "\$SCRATCH_WORK"

# Run ORCA
/cluster/software/commercial/orca/6.0.0/x86_64/bin/orca molecule.inp > molecule.out

# Copy results back to home directory
cp *.gbw *.xyz *.out $HOME_DIR/$basename/

# Cleanup scratch
rm -rf "\$SCRATCH_WORK"
EOF

    # Make script executable
    chmod +x "run_job.sh"
    
    # Submit job
    echo "Submitting job for $basename..."
    sbatch --export=SCRATCH_BASE=$SCRATCH_BASE "run_job.sh"
    
    # Return to base directory
    cd "$HOME_DIR"
}

# Submit jobs for each subdirectory
echo "Submitting calculations for each molecule..."
for dir in "$HOME_DIR"/*/; do
    if [ -d "$dir" ]; then
        molecule=$(basename "$dir")
        echo "Processing $molecule..."
        submit_orca_job "$dir"
    fi
done

echo "All jobs submitted. Check status with 'squeue -u $USER'"
echo "Working directories will be in: $SCRATCH_BASE"
echo "Final output files will be copied back to: $HOME_DIR/[molecule]/"