#!/usr/bin/env python

import sys,getopt
import numpy as np

#def parse_cytof_txt(fp,spt='\t'):
#	f = open(fp)
#	line = f.readline()
#	L = f.readline().strip().split(spt)
#	X = [[float(x) for x in line.strip().split(spt)] for line in f]
#	f.close()
#	return L,np.array(X)

#sometimes data has \r instead of \n
def parse_cytof_txt(fp):
	f = open(fp)
	L = f.read()
	f.close()
	L = [l.split('\t') for l in L.split('\r')]
	X = []
	for l in L[1:]:
		if len(l) == len(L[0]):
			X.append([float(x) for x in l if x])
	return L[0],np.array(X)

def row_normalization(X,p=5):
	y = np.arcsinh(X)
	return y/(np.sum(y**p,axis=-1)**(1./p))[:,np.newaxis]

def write_unfiltered(L,X,fp):
	f = open(fp,'w')
	f.write('\t'.join(['ID']+L)+'\n')
	for i,x in enumerate(X):
		f.write('\t'.join([str(i)]+[str(y) for y in x])+'\n')
	f.close()

def write_filtered(L,X,F,fp,do_norm=False):
	Lf = []
	idx = []
	for i,l in enumerate(L):
		include = True
		for f in F:
			if f in l:
				include = False
				break
		if include:
			Lf.append(l)
			idx.append(i)
	f = open(fp,'w')
	f.write('\t'.join(['ID']+Lf)+'\n')
	if do_norm:
		for i,x in enumerate(row_normalization(X[:,idx])):
			f.write('\t'.join([str(i)]+[str(y) for y in x])+'\n')
	else:
		for i,x in enumerate(X[:,idx]):
			f.write('\t'.join([str(i)]+[str(y) for y in x])+'\n')
	f.close()

if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hi:o:f:',["inputdir=","outputdir=","filters="])
	except:
		print help_message
		sys.exit(2)
	filters = []
	for opt, arg in opts:
		if opt in ('-h','--help'):
			print help_message
			sys.exit()
		elif opt in ('-i','--inputdir'):
			inputdir = arg
		elif opt in ('-o','--outputdir'):
			outputdir = arg
		# filters are ',' separated w/o spaces
		elif opt in ('-f','--filters'):
			filters = arg.split(',')
	L,X = parse_cytof_txt(inputdir)
	write_unfiltered(L,X,outputdir)
	write_filtered(L,X,filters,outputdir+'.filtered')
	write_filtered(L,X,filters,outputdir+'.norm.filtered',do_norm=True)