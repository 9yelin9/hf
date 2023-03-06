#!/bin/bash
. /home/9yelin9/.bash_profile

#$ -q single.q
#$ -pe mpi 1
#$ -j y
#$ -cwd
#$ -o log/$JOB_NAME.log

t0=$(date +%s.%N)
t0_string=$(date)

name="baoso3_dU1.0_half"

for i in `seq 3 1 9`
do
	for eta in `seq 0.1 0.02 0.6`
	do
		bins=$(echo "2^$i"|bc)
		./mod/dos $name p GK $eta $bins
		./mod/dos $name p MK $eta $bins
	done
done

#for eta in `seq 0.02 0.02 0.6`
#do
#	py hf3.py -n $name -ms p GK $eta
#done

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
