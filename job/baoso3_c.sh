#!/bin/bash

#$ -pe mpi 1
#$ -q openmp.q@phase06
#$ -j y
#$ -cwd
#$ -o log/$JOB_NAME.log

source /home/9yelin9/.bash_profile

t0=$(date +%s.%N)
t0_string=$(date)

SOC=0
for ju in `seq 0.0 0.1 0.3`
do
	for n in `seq 0.2 0.2 11.8`
	do
		for u in `seq 0 0.5 5.0`
		do
			./mod/hf3 baoso3 c $ju $SOC $n $u
		done
	done
done

t1=$(date +%s.%N)
t1_string=$(date)

t=$(echo "$t1 - $t0"|bc)
h=$(echo "($t/3600)"|bc)
m=$(echo "($t%3600)/60"|bc)
s=$(echo "($t%3600)%60"|bc)

echo ""
echo "# Job ID       : $JOB_ID"
echo "# Job Name     : $JOB_NAME"
echo "# Time Start   : $t0_string"
echo "# Time End     : $t1_string"
echo "# Time Elapsed : ${h}h ${m}m ${s}s"
echo ""
