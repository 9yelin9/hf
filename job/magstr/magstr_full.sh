#!/bin/bash
. /home/9yelin9/.bash_profile

#$ -q single.q
#$ -pe mpi 1
#$ -j y
#$ -cwd
#$ -o log/$JOB_NAME.log

t0=$(date +%s.%N)
t0_string=$(date)

name="baoso3_dU1.0"

# energy, band, band(gap only)
py hf3.py -n $name -ms e
#py hf3.py -n $name -ms b M
py hf3.py -n $name -ms b G

# ldos
#for eta in `seq 0.1 0.1 0.5`
#do
#	py hf3.py -n $name -ms dl M $eta
#	for i in `seq 3 1 9`
#	do
#		bins=$(echo "2^$i"|bc)
#		./mod/dos $name r ML $eta $bins
#		./mod/dos $name p ML $eta $bins 
#	done
#done

# kdos, peak
for eta in `seq 0.02 0.02 0.4`
do
	#./mod/dos $name d MK $eta
	#./mod/dos $name d MKf $eta
	#./mod/dos $name d MKfb $eta
	./mod/dos $name d GK $eta
	py hf3.py -n $name -ms dd K $eta

	./mod/dos $name p GK $eta 1024
	./mod/dos $name p DK $eta 1024

	for i in `seq 3 1 9`
	do
		bins=$(echo "2^$i"|bc)
		#./mod/dos $name r MK $eta $bins
		#./mod/dos $name r MKf $eta $bins
		#./mod/dos $name r MKfb $eta $bins
		./mod/dos $name r GK $eta $bins
		./mod/dos $name r DK $eta $bins

		#./mod/dos $name p MK $eta $bins
		./mod/dos $name p GK $eta $bins 
		./mod/dos $name p DK $eta $bins
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
