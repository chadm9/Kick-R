#!/bin/tcsh
#PBS -q single
#PBS -l nodes=1:ppn=1
#PBS -l walltime=72:00:00
#PBS -o STDOUT
#PBS -j oe
#PBS -V
#PBS -N Title 
#PBS -A hpc_yexu01
#PBS -m ae
#PBS -M emaillsu.edu

cd /ddnB/work/your_work_directory
python Kick-R.py 
