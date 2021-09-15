#!/bin/sh


#
# Job name
#SBATCH -J dwlp_jf1
#SBATCH -o dex_wt_jf1_%j.out
#SBATCH -e dex_wt_jf1_%j.out

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
export NB_of_jobs=64

# Loop over the job number
frame=50000
arr=(`find . -size 88c |  cut -c 18-25 | paste -sd" "`)


for i in "${arr[@]}"
do
     frame=$((frame+(${i}-1)*5000))
     echo ${frame}
     srun -Q --exclusive --overlap -n 1 -N 1 workScript.sh $i $frame &> worker_${SLURM_JOB_ID}_${i} &
     sleep 1
     frame=50000
done

# keep the wait statement, it is important!

wait
