#!/bin/bash
#SBATCH -t 10:00:00 # walltime,
#SBATCH -N 1 # number of cluster nodes
#SBATCH -n 1 # number of MPI ranks / Slurm tasks
#SBATCH -o slurm_%j.out # pathname of stdout file
#SBATCH -e slurm_%j.err # pathname of stderr file, using job and first node values
#
# load modulefiles
source /usr/projects/climate/SHARED_CLIMATE/anaconda_envs/load_latest_e3sm_unified_grizzly.sh
module purge; module load git; module use /usr/projects/climate/SHARED_CLIMATE/modulefiles/all/
module load gcc/5.3.0 openmpi/1.10.5 netcdf/4.4.1 parallel-netcdf/1.5.0 pio/1.7.2
export MPAS_TOOLS_DIR=/usr/projects/climate/mpeterse/repos/tools/master
export MPAS_SEAICE_TESTCASES_RUN_COMMAND='mpirun -n 1 '
module unload python
source /usr/projects/climate/SHARED_CLIMATE/anaconda_envs/load_latest_e3sm_unified.sh

# run the app
python run_testcase.py
