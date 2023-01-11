#!/bin/bash
. /home/9yelin9/.bash_profile

#$ -q mpi.q
#$ -pe mpi 16
#$ -j y
#$ -cwd
#$ -o log/$JOB_NAME.log

t0=$(date +%s.%N)
t0_string=$(date)

./mod/init cual2o4 f 1 
SAVE="cual2o4"
./mod/init cual2o4 a 1

t1=$(date +%s.%N)
t1_string=$(date)

t=$(echo "$t1 - $t0"|bc)
h=$(echo "($t/3600)"|bc) $SAVE
m=$(echo "($t%3600)/60"|bc)
s=$(echo "($t%3600)%60"|bc)

echo ""
echo "# Job ID       : $JOB_ID"
echo "# Job Name     : $JOB_NAME"
echo "# Time Start   : $t0_string"
echo "# Time End     : $t1_string"
echo "# Time Elapsed : ${h}h ${m}m ${s}s"
echo ""
