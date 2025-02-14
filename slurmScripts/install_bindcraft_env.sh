#!/bin/bash

# package manager
pkg_manager='conda'

install_dir=$(pwd)
SECONDS=0

# check if conda is installed
CONDA_BASE=$(conda info --base 2>/dev/null) || { echo "Error: conda is not installed!!!"; exit 1; }
echo "Conda is installed at: $CONDA_BASE"

# create bindcraft environment
echo "Creating BindCraft environment..."
$pkg_manager create --name BindCraft python=3.10 -y || { echo "Error: Failed to create BindCraft environment"; exit 1; }

# activate environment
echo "Activating the environment..."
source ${CONDA_BASE}/bin/activate ${CONDA_BASE}/envs/BindCraft || { echo "Error: failed to activate BindCraft environment"; exit 1; }
echo "BindCraft environment activated at ${CONDA_BASE}/envs/BindCraft"

# install required packages
echo "Installing required packages..."
$pkg_manager install pip pandas matplotlib numpy"<2.0.0" biopython scipy pdbfixer seaborn libgfortran5 tqdm jupyter ffmpeg pyrosetta fsspec py3dmol chex dm-haiku flax"<0.10.0" dm-tree joblib ml-collections immutabledict optax jaxlib jax cuda-nvcc cudnn -c conda-forge -c nvidia --channel https://conda.graylab.jhu.edu -y || { echo "Error: Failed to install required packages."; exit 1; }

# install colabdesign
pip3 install git+https://github.com/sokrypton/ColabDesign.git --no-deps

# chmod executables
chmod +x ${install_dir}/functions/dssp
chmod +x ${install_dir}/functions/DAlphaBall.gcc

# finish
conda deactivate
printf "BindCraft environment installed\n"

############################################################################################################
############################################################################################################
################## cleanup
printf "Cleaning up ${pkg_manager} temporary files to save space\n"
$pkg_manager clean -a -y
printf "$pkg_manager cleaned up\n"

################## finish script
t=$SECONDS 
printf "Finished setting up BindCraft environment\n"
printf "Activate environment using command: \"$pkg_manager activate BindCraft\""
printf "\n"
printf "Installation took $(($t / 3600)) hours, $((($t / 60) % 60)) minutes and $(($t % 60)) seconds."