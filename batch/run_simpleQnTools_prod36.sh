#!/bin/bash

#SBATCH -D /lustre/stor2/mephi/parfenov/TMP/
#SBATCH -p mephi
#SBATCH --mem-per-cpu=10G
#SBATCH -t 8:00:00
#SBATCH -J QnToolsEfficiencyMpd
#SBATCH -o /lustre/stor2/mephi/parfenov/TMP/slurm_runQnTools_%A_%a.out

list_dir=${1}
output_dir=${2}
ecm=${3}
#nucl_mass=${4}
efficiency_file=${4}
pid_file=${5}

id=$SLURM_ARRAY_TASK_ID

file_list=$( ls $list_dir | head -n $id | tail -n 1 )

mkdir -p $output_dir
cd $output_dir
mkdir $id
cd $id

source /cvmfs/nica.jinr.ru/sw/os/login.sh latest
module add mpddev/v23.09.23-1
export MPDROOT=/lustre/home/user/p/parfenov/Soft/mpdroot/install
source /lustre/home/user/p/parfenov/Soft/mpdroot/install/config/env.sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/home/user/p/parfenov/Soft/QnTools/install-mpd/lib/:/lustre/home/user/p/parfenov/Soft/qntools_macros/build/
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:/lustre/home/user/p/parfenov/Soft/QnTools/install-mpd/include/QnTools/:/lustre/home/user/p/parfenov/Soft/qntools_macros/src/

export exe_correct=/lustre/home/user/p/parfenov/Soft/qntools_macros/build/correct
export exe_correlate=/lustre/home/user/p/parfenov/Soft/qntools_macros/build/correlate

export macro_correct=/lustre/home/user/p/parfenov/Soft/qntools_macros/macro/mpd_fixed_target_correct_prod36.cc
export macro_correlate=/lustre/home/user/p/parfenov/Soft/qntools_macros/macro/mpd_fixed_target_correlate.cc

echo "${exe_correct} ${macro_correct} $list_dir/$file_list $ecm $nucl_mass $efficiency_file $pid_file"

# PLAIN
time $exe_correct $macro_correct $list_dir/$file_list $ecm $efficiency_file $pid_file
# RECENTERING
time $exe_correct $macro_correct $list_dir/$file_list $ecm $efficiency_file $pid_file
# TWIST AND RESCALING
time $exe_correct $macro_correct $list_dir/$file_list $ecm $efficiency_file $pid_file

echo "${exe_correlate} ${macro_correlate} correction_out.root"
time ${exe_correlate} ${macro_correlate} correction_out.root

echo "The End."
