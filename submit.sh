#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=112
#SBATCH --time=01:30:00
#SBATCH -J compile_WRF
#SBATCH --account=TG-ATM170028
#SBATCH --partition=spr
#SBATCH --mail-type=ALL
#SBATCH --mail-user=storm@mit.edu

source "source_this_file.sh"
wait
echo $?
#./clean -a
#wait
#echo $?
echo -e "78\n1" | ./configure
wait
echo $?
./compile em_les >& log.em_les
