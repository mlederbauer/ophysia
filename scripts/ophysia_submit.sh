#!/bin/bash

# Check if directory name is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 DIRNAME"
    echo "Example: $0 my_calculation"
    exit 1
fi

# Get directory name from command line argument
CALC_DIR="$1"

# Directory setup
HOME_DIR="/cluster/home/$USER/$CALC_DIR"
SCRATCH_DIR="/cluster/scratch/$USER/$CALC_DIR"

# Default job parameters
MEM_PER_CPU="4000"
NTASKS="2"

# Base directories
HOME_DIR="$PWD/$CALC_DIR"
SCRATCH_BASE="/cluster/scratch/$USER/$CALC_DIR"

# Ensure scratch base directory exists
mkdir -p "$SCRATCH_BASE"

# Load required modules
module purge
module load stack/2024-06
module load gcc/12.2.0
module load orca/6.0.0
module load openmpi/4.1.6

# Function to submit a job
submit_orca_job() {
    local dir=$1
    local basename=$(basename "$dir")
    local scratch_dir="$SCRATCH_BASE/$basename"
    
    # Set walltime based on directory name
    if [[ $basename == *go_sp ]]; then
        WALLTIME="24:00:00"
    else
        WALLTIME="4:00:00"  # Default walltime
    fi
    
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
cp *.gbw *.xyz *.out *.hess $HOME_DIR/$basename/

# Cleanup scratch
rm -rf "\$SCRATCH_WORK"
EOF

    # Make script executable
    chmod +x "run_job.sh"
    
    # Submit job
    echo "Submitting job for $basename (walltime: $WALLTIME)..."
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