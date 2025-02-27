#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --partition=gen_gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:Tesla_V100-PCIE-32GB:1
#SBATCH --mem 42gb
#SBATCH --time 24:00:00
#SBATCH --job-name=run_juptNote
#SBATCH --output=/hpf/projects/mtyers/ningrui/NXBindCraft/Logs/jupLogs/jupyter_%j.out


# Load module (not necessary if using local bindcraft env)
#echo "Loading BindCraft module..."
#module load BindCraft/

# Initialise environment and modules
CONDA_BASE=/hpf/projects/mtyers/ningrui/miniconda3
source ${CONDA_BASE}/bin/activate BindCraft # ${CONDA_BASE}/envs/BindCraft
export LD_LIBRARY_PATH=${CONDA_BASE}/lib

echo "Python in use:"
which python

# Start the Jupyter Notebook server
jupyter notebook --no-browser --ip=0.0.0.0 --port=8878