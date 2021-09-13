#!/bin/sh
#
# Job name
#SBATCH -J dw_jf
#SBATCH -o dex_wt_jf_%j.out
#SBATCH -e dex_wt_jf_%j.out

##SBATCH --mail-user=yulian.gavrilov@bpc.lu.se
##SBATCH --mail-type=ALL

#SBATCH -A LU2021-2-64
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=3100
#SBATCH -N 1
#SBATCH --cpus-per-task=5
#
# Time needed to complete the job
#SBATCH -t 168:00:00
#

cat $0

# set the number of jobs - change for your requirements
export NB_of_jobs=5

# Loop over the job number
frame=50000

for ((i=1; i<=$NB_of_jobs; i++))
do
    srun -Q --exclusive --overlap -n 1 -N 1 workScript.sh $i $frame &> worker_${SLURM_JOB_ID}_${i} &
    sleep 1
    frame=$((frame+5000))
done

# keep the wait statement, it is important!

wait
