#!/bin/bash

#SBATCH --job-name gadma_Larix
#SBATCH -A inbreh
#SBATCH -t 7-00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=53
#SBATCH --mem=300G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/inbreh/turck_fsc/errs_outs/err_Larix_gadma_%A_%a.err
#SBATCH -o /project/inbreh/turck_fsc/errs_outs/std__Larix_gadma_%A_%a.out


cd /project/inbreh/turck_fsc/gadma   # set the correct working directory 

# load modules and conda environment
module load miniconda3/24.3.0
conda activate gadma

# Run gadma
gadma -p gadma_Larix_K3_2event_params.txt 
