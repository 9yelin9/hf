#!/bin/bash
. /home/9yelin9/.bash_profile

#$ -q single.q
#$ -pe mpi 1
#$ -j y
#$ -cwd
#$ -o log/$JOB_NAME.log

t0=$(date +%s.%N)
t0_string=$(date)

py hf3.py -n baoso3 -ms d 1m-0-0 100 0.5 1.0
py hf3.py -n baoso3 -ms d 1f-0-0 100 0.5 1.0
py hf3.py -n baoso3 -ms d 1b-0-0 100 0.5 1.0
py hf3.py -n baoso3 -ms d 1fb-0-0 100 0.5 1.0

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
