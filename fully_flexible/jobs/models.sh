#!/usr/bin/bash

#SBATCH -c 1 ## number of cores
#SBATCH -t 0-03:00 ## amount of time in D-HH:MM
#SBATCH -p shared ## Partition to submit to 
#SBATCH --mem=8000 ## memory pool for all cores
#SBATCH -o logs/models/log.stdout_%a ## STDOUT 
#SBATCH -e logs/models/log.stderr_%a ## STDERR
#SBATCH --account=haneuse_lab
#SBATCH --array=1-1000

module load R/4.2.2-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER

cd $HOME/Levis_Phd_Paper1-ExtraSims/fully_flexible
  
Rscript estimate_nuisance_fx.R $SLURM_ARRAY_TASK_ID
