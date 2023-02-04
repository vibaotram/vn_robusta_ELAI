#!/bin/bash
#SBATCH --job-name group_fst
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --partition=highmem
#SBATCH --error /data3/projects/vietcaf/baotram/scripts/robusta_vn/test_10_clones/test_elai/genetic_structure/slurm-%x_%j.log
#SBATCH --output /data3/projects/vietcaf/baotram/scripts/robusta_vn/test_10_clones/test_elai/genetic_structure/slurm-%x_%j.log

module load system/singularity/3.6.0
singularity exec /home/baotram/singularity-container_myr_4-0-2_rstudio_1.3.sif Rscript /data3/projects/vietcaf/baotram/scripts/robusta_vn/test_10_clones/test_elai/genetic_structure/group_fst.R