#!/bin/bash
#SBATCH --job-name="drug_simulation"
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=6
#SBATCH --cpus-per-task=8
#SBATCH -t 22:00:00

/apps/GREASY/latest/INTEL/IMPI/bin/greasy run_drug_simulations.sh
