#!/bin/bash
#SBATCH --job-name=modules_job  # A short name for your job
#SBATCH --nodes=1               # Request one node
#SBATCH --ntasks=1              # Run a single task
#SBATCH --cpus-per-task=1       # Number of CPU cores per task
#SBATCH --mem=64G               # **Requested memory (e.g., 64 Gigabytes)**
#SBATCH --time=04:00:00         # Max job run time (HH:MM:SS)
#SBATCH --output=slurm-%j.out   # Standard output log (optional)
#SBATCH --error=slurm-%j.err    # Standard error log (optional)
#SBATCH --mail-type=END,FAIL    # Email alerts when job ends or fails
#SBATCH --mail-user=ym3292@princeton.edu

# Load the specific Anaconda module
module load anaconda3/2024.10

# Activate your Conda environment
source activate qcb_env_adroit

# Navigate to your script's location (optional if already there)
# cd /scratch/network/ym3292/coessentiality

# Execute your Python script
python PCA_visualize.py

# Deactivate the environment (optional)
# conda deactivate