import os
import re
import argparse
import numpy as np

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-st', '--strain', help='<strain>')
parser.add_argument('-du',  '--dU',  type=float, help='<dU>')
parser.add_argument('-uf',  '--UF',  type=int,   help='<UF>')

parser.add_argument('--band' type=int,   help='<Nkb>')
parser.add_argument('--DOS', type=float, help='<ep>')
args = parser.parse_args()                                                                     

fd = open('job/default.txt', 'r')
save = '%s_dU%.1f_UF%d' % (args.strain, args.dU, args.UF)
os.makedirs('output/%s' % save, exist_ok=True)

if   args.BAND: save_job = save + '_band_Nkb%d' % args.BAND
elif args.DOS:  save_job = save + '_dos_ep%.2f' % args.DOS
else:           save_job = save

os.makedirs('job/%s' % save_job, exist_ok=True)
os.makedirs('log/%s' % save_job, exist_ok=True)

js = np.arange(0, 0.3, 0.1)
ts = ['F0'] + ['A0', 'A2', 'A5'] + ['C0', 'C1', 'C6'] + ['G0']

for t in ts:
	for j in js:
		fn = 'job/%s/%s_JU%.2f.sh' % (save_job, t, j)
		f = open(fn, 'w')

		for line in fd:
			if re.search('JOB_NAME[.]log', line):
				f.write('#$ -o log/%s/$JOB_NAME.log\n' % save_job)

			elif re.search('###', line):
				if re.search('F', t): f.write('for n in `seq 0.1 0.1 5.9`\ndo\n')
				else:                 f.write('for n in `seq 0.2 0.2 11.8`\ndo\n')

				f.write('\tfor u in `seq 0 %.1f %d`\n\tdo\n' % (args.dU, args.UF))

				if   args.BAND: f.write('\t\t./hf/hf %s %s %.2f 0 $n $u %d\n\tdone\ndone\n'   % (save, t, j, args.BAND))
				elif args.DOS:  f.write('\t\t./hf/hf %s %s %.2f 0 $n $u %.2f\n\tdone\ndone\n' % (save, t, j, args.DOS))
				else:           f.write('\t\t./hf/hf %s %s %.2f 0 $n $u\n\tdone\ndone\n'      % (save, t, j))
			else: f.write(line)

		print(fn)
		f.close()
		fd.seek(0)
fd.close()
