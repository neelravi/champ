#!/bin/bash -x 
#SBATCH --account=prcoe10
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -J c60input
#SBATCH -o c60.%j.out
#SBATCH -e c60.%j.err
#SBATCH --exclusive
#SBATCH --partition=batch

##################################################################################
# Load the modules to create the proper environment
##################################################################################
###source /p/home/jusers/landinez1/juwels/Codes/pyscf_env/bin/activate
###module load Stages/2023
###module load GCC/11.3.0
###module load OpenMPI/4.1.4
###module load ScaLAPACK/2.2.0-fb


###module load Stages/2023
###module load Intel/2022.1.0
###module load IntelMPI/2021.6.0
###module load imkl/2022.1.0 

module load Stages/2023
module load GCC/11.3.0
module load PySCF/2.1.1


##################################################################################

cd $PWD

srun --cpus-per-task 48 --cpu-bind=cores --ntasks 1 python C60.py
sleep 1s


