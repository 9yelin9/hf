#!/bin/bash
. /home/9yelin9/.bash_profile

#$ -q openmp.q@phase09
#$ -pe mpi 1
#$ -j y
#$ -cwd
#$ -o log/baoso3_dU1.0/$JOB_NAME.log

t0=$(date +%s.%N)
t0_string=$(date)

for n in `seq 0.2 0.2 11.8`
do
	for u in `seq 0 1.0 8`
	do
		./mod/hf3 baoso3 baoso3_dU1.0 0.00 0 c1 $n $u
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
