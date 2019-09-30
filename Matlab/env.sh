#! /usr/bin/bash

if [[ "${BASH_SOURCE[0]}" != "${0}" ]]
then
	echo "Script is being sourced"
else
	echo "ERROR: Script is a subshell"
	echo "To affect your current shell enviroment source this script with:"
	echo "source env.sh"
	exit
fi

# check that this file is run
export ENV_RUN=true

# if hpipm folder not specified assume alongside the parent of this folder
HPIPM_MAIN_FOLDER=${HPIPM_MAIN_FOLDER:-"$(pwd)/../../hpipm"}
export HPIPM_MAIN_FOLDER
echo
echo "HPIPM_MAIN_FOLDER=$HPIPM_MAIN_FOLDER"

# if blasfeo folder not specified assume alongside the parent of this folder
BLASFEO_MAIN_FOLDER=${BLASFEO_MAIN_FOLDER:-"$(pwd)/../../blasfeo"}
export BLASFEO_MAIN_FOLDER
echo
echo "BLASFEO_MAIN_FOLDER=$BLASFEO_MAIN_FOLDER"

# export matlab_octave_mex folder
# MATLAB case
export MATLABPATH=$MATLABPATH:$HPIPM_MAIN_FOLDER/interfaces/matlab_octave/
echo
echo "MATLABPATH=$MATLABPATH"
# Octave case
export OCTAVE_PATH=$OCTAVE_PATH:$HPIPM_MAIN_FOLDER/interfaces/matlab_octave/
echo
echo "OCTAVE_PATH=$OCTAVE_PATH"

#export HPIPM_MEX_FLAGS="GCC=/usr/bin/gcc-4.9"

# export LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HPIPM_MAIN_FOLDER/lib:$BLASFEO_MAIN_FOLDER/lib
echo
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

