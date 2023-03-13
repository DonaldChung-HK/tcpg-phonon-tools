#! /bin/bash -l
#SBATCH -p scarf
#SBATCH --job-name 25_diiodothiophene_CASTEP_opt
#SBATCH --nodes=6
#SBATCH --exclusive
#SBATCH -C amd
#SBATCH --time=10:00:00
#SBATCH --array=0 # job array index remove for singleJob

module purge
module load AMDmodules
module load Python/3.10.4-GCCcore-11.3.0
module load CASTEP/21.1.1-iomkl-2021a
module list

mpirun castep.mpi CuSO4_5H2O

#find . -name "*.check" -delete