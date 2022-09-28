#!/bin/bash

#$ -pe mpi 16
#$ -q openmp.q
#$ -j y
#$ -cwd

export PATH="/home/9yelin9/usr/bin:$PATH"
export LD_LIBRARY_PATH="/usr/lib:/usr/lib64:/usr/slib:/home/9yelin9/usr/lib:/home/9yelin9/usr/lib64:$LD_LIBRARY_PATH"

module load openmpi/gcc-4.8.5/4.1.0
module load gsl/gcc-4.8.5/2.7.1

t0=$(date +%s.%N)
t0_string=$(date)

SOC=0
for ju in `seq 0.1 0.1 0.3`
do
	for n in `seq 0.2 0.2 11.8`
	do
		for u in `seq 0 1 9`
		do
			./baoso3 g $ju $SOC $n $u
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
echo "# Job Start : $t0_string"
echo "# Job End   : $t1_string"
echo "# Elapsed time : ${h}h ${m}m ${s}s"
