#!/bin/bash

#SBATCH -D /lustre/stor2/mephi/parfenov/TMP/
#SBATCH -p mephi
#SBATCH --mem-per-cpu=10G
#SBATCH -t 8:00:00
#SBATCH -J QnToolsNeutronBmn
#SBATCH -o /lustre/stor2/mephi/parfenov/TMP/slurm_runQnToolsBmn_%A_%a.out
#SBATCH -e /lustre/stor2/mephi/parfenov/TMP/slurm_runQnToolsBmn_%A_%a.out

list_dir=${1}
output_dir=${2}

id=$SLURM_ARRAY_TASK_ID

file_list=$( ls $list_dir | head -n $id | tail -n 1 )

mkdir -p $output_dir
cd $output_dir
mkdir $id
cd $id

source /cvmfs/bmn.jinr.ru/config/x86_64-centos7/cluster_config.sh
source /lustre/home/user/p/parfenov/Soft/bmnroot_jul2024/build/config.sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/home/user/p/parfenov/Soft/QnTools/install-bmn/lib/:/lustre/home/user/p/parfenov/Soft/qntools_macros/build_bmn/
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:/lustre/home/user/p/parfenov/Soft/QnTools/install-bmn/include/QnTools/:/lustre/home/user/p/parfenov/Soft/qntools_macros/src/

export exe_correct=/lustre/home/user/p/parfenov/Soft/qntools_macros/build_bmn/correct
export exe_correlate=/lustre/home/user/p/parfenov/Soft/qntools_macros/build_bmn/correlate

export macro_correct=/lustre/home/user/p/parfenov/Soft/qntools_macros/macro/jam_neutron_correct.cc
export macro_correlate=/lustre/home/user/p/parfenov/Soft/qntools_macros/macro/jam_neutron_correlate.cc

echo "${exe_correct} ${macro_correct} $list_dir/$file_list"

# PLAIN
time $exe_correct $macro_correct $list_dir/$file_list
# RECENTERING
time $exe_correct $macro_correct $list_dir/$file_list
# TWIST AND RESCALING
time $exe_correct $macro_correct $list_dir/$file_list

echo "${exe_correlate} ${macro_correlate} correction_out.root"
time ${exe_correlate} ${macro_correlate} correction_out.root

echo "The End."
