#!/bin/bash
#SBATCH --job-name validate_elai
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --partition=long
#SBATCH --error /data3/projects/vietcaf/baotram/scripts/robusta_vn/test_10_clones/test_elai/validate_elai/slurm-%x_%j.log
#SBATCH --output /data3/projects/vietcaf/baotram/scripts/robusta_vn/test_10_clones/test_elai/validate_elai/slurm-%x_%j.log

module load system/singularity/3.6.0
singularity exec -B /home/baotram /home/baotram/singularity-container_myr_4-0-2_rstudio_1.3.sif \
Rscript /data3/projects/vietcaf/baotram/scripts/robusta_vn/test_10_clones/test_elai/validate_elai/validate_elai.R