#!/bin/tcsh
#PBS -q single
#PBS -l nodes=1:ppn=2
#PBS -l walltime=72:00:00
#PBS -o STDOUT
#PBS -j oe
#PBS -V
#PBS -N Title 
#PBS -A hpc_yexu01 
#PBS -m ae
#PBS -M email@lsu.edu

set NPROCS=`wc -l $PBS_NODEFILE |gawk '//{print $1}'`
set LINDA=`cat $PBS_NODEFILE | uniq | tr '\\n' "p" |sed 's/p/:2,/g' | sed 's|,$||' `
setenv GAUSS_SCRDIR /ddnB/work/your_scratch_directory
source $g09root/g09/bsd/g09.login
cd /ddnB/work/your_work_directory 

foreach i (*.kick)
cat $i| sed "s/LINDA/$LINDA/" > temp$$.inp
g09 < temp$$.inp > $i.out
rm -f temp$$.inp
end

