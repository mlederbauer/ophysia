#!/bin/bash

# Define job IDs directly in the script
job_ids=(
    "ethanol_sp_only"
    "ethanol_go_sp"
    "gabaarg_go_sp"
    "gbaarg_sp_only"
    "id_3_go_sp"
    "id_3_sp_only"
    "id_108_go_sp"
    "id_108_sp_only"
    "id_123_go_sp"
    "id_123_sp_only"
)

# Check if ophysia_submit.sh exists and is executable
if [ ! -x "./ophysia_submit.sh" ]; then
    echo "Error: ophysia_submit.sh not found or not executable in current directory"
    exit 1
fi

# Execute jobs one by one
for jobid in "${job_ids[@]}"; do
    echo "Submitting job: $jobid"
    ./ophysia_submit.sh "$jobid"
    
    # Add a small delay between submissions (optional)
    sleep 1
done

echo "All jobs submitted successfully"