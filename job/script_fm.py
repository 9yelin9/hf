import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import re
import argparse
import numpy as np

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-du', '--dU', type=float, help='<dU>')
args = parser.parse_args()                                                                     

fd = open('job/default.txt', 'r')
save = 'FM_dU%.1f' % args.dU
os.makedirs('job/%s'    % save, exist_ok=True)
os.makedirs('log/%s'    % save, exist_ok=True)
os.makedirs('output/%s' % save, exist_ok=True)

js = np.arange(0, 0.3, 0.1)
ts = ['F0']

os.makedirs('job/%s' % save, exist_ok=True)
os.makedirs('log/%s' % save, exist_ok=True)
os.makedirs('output/%s' % save, exist_ok=True)

for t in ts:
	for j in js:
		f = open('job/%s/%s_JU%.2f.sh' % (save, t, j), 'w')

		for line in fd:
			if re.search('JOB_NAME[.]log', line):
				f.write('#$ -o log/%s/$JOB_NAME.log\n' % save)

			elif re.search('###', line):
				f.write('for n in `seq 0.1 0.1 5.9`\ndo\n')
				f.write('\tfor u in `seq 0 %.1f 8`\n\tdo\n' % args.dU)
				f.write('\t\t./main/hf %s %s %.2f 0 $n $u\n\tdone\ndone\n' % (save, t, j))
			else: f.write(line)

		f.close()
		fd.seek(0)
fd.close()
