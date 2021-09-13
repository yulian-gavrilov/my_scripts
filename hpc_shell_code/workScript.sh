#!/bin/sh
# document this script to stdout (assumes redirection from caller)
cat $0

# receive my worker number
export WRK_NB=$1
export WRK_frame=$2

# create worker-private subdirectory in tmp
export EXE_DIR=$PWD
export WRK_DIR=$PWD/WRK_${WRK_NB}
mkdir $WRK_DIR

# create a variable to address the "job directory"
export JOB_DIR=$PWD/md5ns_runs/job_${WRK_NB}

# now copy the input data and program from there

cd $JOB_DIR

cp -p * $WRK_DIR

# change to the execution directory

cd $WRK_DIR

# run the program

module add GCC/8.3.0  CUDA/10.1.243  OpenMPI/3.1.4
module add GROMACS/2019.4-PLUMED-2.5.4

gmx mdrun  -nt 5 -ntomp 5 -deffnm md5ns_frame${WRK_frame} -v -cpi

# rescue the results back to job directory

cp -p * ${JOB_DIR}

# clean up the local disk and remove the worker-private directory

cd $EXE_DIR

rm -rf WRK_${WRK_NB}
