#!/bin/bash -l
#SBATCH --job-name=gwas_gr20_umax
#SBATCH --account=project_2000350
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ilkka.kronholm@jyu.fi

# Load r-env-singularity
module load r-env-singularity

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2000350" >> ~/.Renviron

# Match thread and core numbers
export SINGULARITYENV_OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Thread affinity control
export SINGULARITYENV_OMP_PLACES=cores
export SINGULARITYENV_OMP_PROC_BIND=close

srun singularity_wrapper exec Rscript --no-save asso_gr20_umax.R
