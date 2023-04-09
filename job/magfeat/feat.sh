#!/bin/bash
. /home/9yelin9/.bash_profile

#$ -q openmp.q@phase06
#$ -pe mpi 1
#$ -j y
#$ -cwd
#$ -o log/magfeat/$JOB_NAME.log

t0=$(date +%s.%N)
t0_string=$(date)

for eta in `seq 0.1 0.1 0.3`
do
	./mod/feat baoso3_dU1.0 d MK $eta
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
