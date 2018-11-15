#!/bin/bash
#PBS -N gpuMATLAB_afh115
#PBS -l walltime=00:35:00
#PBS -l select=1:ncpus=8:mem=16gb:ngpus=1:gpu_type=K80
												# K80 for more power, P1000 for chill times
												# Can check availability using 'availability' command in bash terminal

# # Clone KiloSort and npy-matlab in your code folder
# git clone https://github.com/cortex-lab/KiloSort.git
# git clone https://github.com/kwikteam/npy-matlab.git

srcDir=$HOME # Folder where your data is / master_XYZ.m is
# pwd # States current working directory
# echo CUDA Visible Devices: $CUDA_VISIBLE_DEVICES  End List

# Load New GCC
	module load gcc/6.2.0 # Needed for compatibility with CUDA 9
# Load CUDA
	module load cuda/9.0.176
# Load MATLAB/VERSION
	module load matlab/R2018a
# Load KiloSort and numpy required
	cp -R $HOME/KiloSort/* $TMPDIR/
# Make mex use C++
	mex -setup C++ -v # Verbose just so you can make sure it all good
#Include the .mexa64 files
	PATH=$PATH:/$TMPDIR/KiloSort/CUDA/
	export PATH=$PATH
	echo $PATH
	# ls $TMPDIR/KiloSort/CUDA/
# Compile the mexGPU code
	# matlab -nodisplay -nodesktop -r "run $TMPDIR/KiloSort/CUDA/mexGPUall.m , exit"
		# Platform dependent, as we are using Linux, compiles .mexa64 files
		# Should work on all cards, if not just uncomment this line (will then compile each time it runs)
			# Might also neeed to remove the existing .mexa64 files in the folder from above
# 

# List TMPDIR Structure (to help debugging)
# ls $TMPDIR -Ra

FILEPATH=$TMPDIR/KiloSort/eMouse 	# Assume you are already in the $TMPDIR
FILENAME="master_eMouse"
matlab -nodisplay -nodesktop -r "run $FILEPATH/$FILENAME.m, exit"

# Storing all processed files into directory
# $WORK a WORK directory dedicated to the user: /work/afh115
	mkdir $WORK/$PBS_JOBID
	cp -a * $WORK/$PBS_JOBID

# Stores location of output directory to file, '>>' means appends to end of file
  date >> $HOME/outputLocationGPUjob.txt 
  echo $TMPDIR/ >> $HOME/outputLocationGPUjob.txt 
  
# Copies processed data to srcDir
  mkdir $srcDir/eMouseTestDataOutput
	cp -a * $srcDir/eMouseTestDataOutput/
 # Removes the extra stuff we copied
  rm -rf KiloSort/
  rm -rf npy-matlab/