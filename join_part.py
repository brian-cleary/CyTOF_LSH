#!/usr/bin/env python

import sys,getopt
import glob,os
from collections import defaultdict
from operator import itemgetter


#header="(mclheader\nmcltype matrix\ndimensions NUMCELLSxNUMCELLS\n)\n(mclmatrix\nbegin\n"

if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'i:o:r:t:k:',["inputdir=","outputdir=","rank=","thresh=","knearest="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-i','--inputdir'):
			inputdir = arg
		elif opt in ('-o','--outputdir'):
			outputdir = arg
		elif opt in ('-r','--rank'):
			fr = int(arg) -1
		elif opt in ('-t','--thresh'):
			thresh = float(arg)
		elif opt in ('-k','--knearest'):
			knearest = int(arg)
	FP = glob.glob(os.path.join(inputdir+'/distances','*.distances.tab'))
	D = defaultdict(dict)
	x_min = fr*1000
	x_max = (fr + 1)*1000
	largest_n = 0
	for fp in FP:
		f = open(fp)
		for line in f:
			ls = line.strip().split('\t')
			n = int(ls[0])
			if n >= x_min:
				if n < x_max:
					for x in ls[1:]:
						try:
							i = x.split(':')
							largest_n = max(int(i[0]),largest_n)
							if float(i[1]) > thresh:
								D[n][int(i[0])] = float(i[1])
						except:
							pass
				else:
					break
		f.close()
	if len(D) > 0:
		f = open(inputdir+'sampled_data.tab.norm.filtered')
		num_cols = len(f.readline().split('\t'))-1
		f.close()
		f = open('%s/all_data.%d' % (outputdir,fr),'w')
		K = D.keys()
		K.sort()
		# knearest edges per 10000 cells
		# this will retain more edges for larger graphs
		# BADNESS (goodness?) WARNING: grows linearly with the number of cells
		knearest = knearest*largest_n/10000
		for k in K:
			d = []
			for a,b in sorted(D[k].iteritems(),key=itemgetter(1),reverse=True):
				if b < 1:
					if len(d) < knearest:
						d.append((a,b))
					else:
						break
				else:
					d.append((a,b))
			d.sort()
			d = [str(a)+':'+str(b) for a,b in d]
			f.write('%d:%d '% (k,num_cols) + ' '.join(d)+' $\n')
			del D[k]
		#f.write(')\n')
		f.close()