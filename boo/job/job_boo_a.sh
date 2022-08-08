#!/bin/bash

#$ -pe mpi 1
#$ -q openmp.q
#$ -j y
#$ -cwd

export PATH="/home/9yelin9/usr/bin:/home/9yelin9/.local/bin:$PATH"
export LD_LIBRARY_PATH="/usr/lib:/usr/lib64:/usr/slib:/home/9yelin9/usr/lib:/home/9yelin9/usr/lib64:/home/9yelin9/.local/lib:$LD_LIBRARY_PATH"

t0=$(date +%s.%N)
t0_string=$(date)

JU=0
SOC=0
IS_UNFOLD=0
for n in `seq 0.2 0.2 11.8`
do
	for u in `seq 0 1 9`
	do
		./boo a $JU $SOC $n $u $IS_UNFOLD
	done
done

t1=$(date +%s.%N)
t1_string=$(date)

t=$(echo "$t1 - $t0"|bc)
h=$(echo "($t/3600)"|bc)
m=$(echo "($t%3600)/60"|bc)
s=$(echo "($t%3600)%60"|bc)

echo ""
echo "# Job Start : $t0_string"
echo "# Job End   : $t1_string"
echo "# Elapsed time : ${h}h ${m}m ${s}s"
