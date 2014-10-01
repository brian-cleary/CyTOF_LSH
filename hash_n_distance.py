#!/usr/bin/env python

import sys,getopt,os
import numpy as np
import random
from scipy.spatial import distance

def hyperplane_hash(fp,depth=5):
	if '.npy' in fp:
		X = np.load(fp)
	else:
		f = open(fp)
		X = []
		line = f.readline()
		for line in f:
			X.append([float(x) for x in line.strip().split('\t')[1:]])
		f.close()
		X = np.array(X)
	HashedValues = np.ones(len(X),dtype=np.int64)
	HashedKeys = []
	RemainingHashes = {1: True}
	while len(RemainingHashes) > 0:
		H = RemainingHashes.keys()
		for h in H:
			ix = np.where(HashedValues == h)[0]
			if len(ix) > depth*X.shape[1]:
				x = X[ix]
				p = one_hyperplane(x)
				L = (np.dot(x,p[:-1]) > p[-1])
				HashedValues[ix] += h + L
				RemainingHashes[2*h] = True
				RemainingHashes[2*h + 1] = True
			HashedKeys.append(h)
			del RemainingHashes[h]
			#print h,len(ix)
	return X,HashedValues,HashedKeys

def one_hyperplane(x):
	y = x[random.sample(xrange(x.shape[0]),x.shape[1])]
	z = np.zeros((y.shape[0]+1,y.shape[1]+1))
	z[:-1,:-1] = y
	z[:,-1] = -1
	U,W,V = np.linalg.svd(z)
	return V[-1,:]

if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'i:o:r:h:p:d:',["inputdir=","outputdir=","rank=","hashes=","depth="])
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
		elif opt in ('-h','--hashes'):
			h = int(arg)
		elif opt in ('-p'):
			p = int(arg)
		elif opt in ('-d','--depth'):
			if '.' in arg:
				depth = float(arg)
			else:
				depth = int(arg)
	X,HV,HK = hyperplane_hash(inputdir,depth=depth)
	outfile = '%s/%d.distances.tab' % (outputdir,fr)
	f = open(outfile,'w')
	for i in HK:
		ix = np.where(HV == i)[0]
		if len(ix) > 0:
			x = X[ix]
			D = 1 - distance.squareform(distance.pdist(x,'minkowski',p))
			for j,d in enumerate(D):
				f.write('\t'.join([str(ix[j])] + [str(ix[k])+':'+str(d[k]) for k in range(len(d)) if d[k]>0.25]) + '\n')
	f.close()
	del X
	os.system('sort -nk1,1 %s -o %s' % (outfile,outfile))
		