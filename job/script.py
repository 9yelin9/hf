import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import re
import argparse
import numpy as np

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
#parser.add_argument('-n', '--name', type=str, default=None, help='<name of save dir>')
parser.add_argument('-q', '--queue', type=int, default=9, help='<num of queue>')
args = parser.parse_args()                                                                     

fd = open('job/default.txt', 'r')

#save = '%s' % (args.name)
name = 'baoso3_dU0.2'
save = 'spec'

os.makedirs('job/%s' % save, exist_ok=True)
os.makedirs('log/%s' % save, exist_ok=True)

for eta in np.arange(0.1, 0.31, 0.05):
	f = open('job/%s/%s_eta%.2f.sh' % (save, save, eta), 'w')
	for line in fd:
		if re.search('openmp[.]q', line):
			f.write('#$ -q openmp.q@phase%02d\n' % args.queue)
		elif re.search('JOB_NAME[.]log', line):
			f.write('#$ -o log/%s/$JOB_NAME.log\n' % save)
		elif re.search('###', line):
			f.write('./mod/spec %s s M %.2f\n' % (name, eta))
		else: f.write(line)
	f.close()
	fd.seek(0)
fd.close()

"""
dU = float(re.sub('dU', '', re.search('dU\d+[.]\d+', save).group()))

#js = [0.5*i for i in range(5)]
js = [0.1*i for i in range(3)]
#js = [0.1*i + 0.05 for i in range(2)]
ts = ['a0', 'a2', 'a5'] + ['c0', 'c1', 'c6'] + ['g0']

os.makedirs('job/%s' % save, exist_ok=True)
os.makedirs('log/%s' % save, exist_ok=True)
os.makedirs('output/%s' % save, exist_ok=True)

for j in js:
	for t in ts:
		f = open('job/%s/%s_JU%.2f_%s.sh' % (save, save, j, t), 'w')

		for line in fd:
			if re.search('openmp[.]q', line):
				f.write('#$ -q openmp.q@phase%02d\n' % args.queue)

			elif re.search('JOB_NAME[.]log', line):
				f.write('#$ -o log/%s/$JOB_NAME.log\n' % save)

			elif re.search('###', line):
				f.write('for n in `seq 0.2 0.2 11.8`\ndo\n')
				#f.write('for n in "6.0"\ndo\n')
				f.write('\tfor u in `seq 0 %.1f 8`\n\tdo\n' % dU)
				#f.write('\tfor u in `seq 0 0.2 8`\n\tdo\n')
				f.write('\t\t./mod/hf3 baoso3 %s %.2f 0 %s $n $u\n\tdone\ndone\n' % (save, j, t))
			else: f.write(line)

		f.close()
		fd.seek(0)
fd.close()
"""
