#!/bin/bash
#
# Installation script example on the LUCIA Tier-1 facility of the cenaero
#
# Tested with the following modules 
# module load EasyBuild/2022a
# module load OpenMPI/4.1.4-GCC-11.3.0
# module load FFTW/3.3.10-GCC-11.3.0
# module load OpenBLAS/0.3.20-GCC-11.3.0
#
# We follow the instruction from the CECI to install custom softwares from sources in a ~/.local directory
# https://support.ceci-hpc.be/doc/_contents/UsingSoftwareAndLibraries/CompilingSoftwareFromSources/index.html
# This is optional if you already have custom software installed. Do not forget to add a path to the ~/.local directory (last step)
#
# BeamDynLib is then compiled with the Makefile.lucia file and installed. 


# Error handling (stops the scripts if something goes wrong)
handle_error() {
    echo "An error occurred on line $1"; echo "Installation FAILED :("
    exit 1
}

trap 'handle_error $LINENO' ERR



# Create directory for custom software installed from sources
cd $HOME
mkdir -p ~/.local/{bin,lib,src,include}
cd ~/.local/src

# Clone, compile and install BeamDynLib
git clone https://git.immc.ucl.ac.be/BigFlow/beamdyn-library BeamDynLib
cd BeamDynLib
mv Makefile Makefile.darwin
mv Makefile.lucia Makefile
make -j 8
make install

# Add the newly install software to the environment  (! This is a manual step !)
echo
echo
echo 'Please add the following to your ~/.bashrc file'
echo 'export PATH=~/.local/bin:$PATH'
echo 'export LD_LIBRARY_PATH=~/.local/lib:$LD_LIBRARY_PATH'
echo 'export LIBRARY_PATH=~/.local/lib:$LIBRARY_PATH'
echo 'export CPATH=~/.local/include:$CPATH'
echo
echo 'Installation SUCCESS :)'
