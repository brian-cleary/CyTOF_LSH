#!/usr/bin/env python

import sys, getopt, os, glob
import numpy as np
from collections import defaultdict
import Pycluster
from scipy import misc
from operator import itemgetter
from scipy.spatial import distance
from multiprocessing import Pool

def get_data(fp):
	f = open(fp)
	header = f.readline().strip().split('\t')
	Channels = [s.split(': ')[1] for s in header[1:]]
	Ids = {}
	X = []
	for line in f:
		ls = line.strip().split('\t')
		Ids[int(ls[0])] = len(X)
		X.append([float(x) for x in ls[1:]])
	f.close()
	return Channels,Ids,np.array(X)

def get_alpha_thresholds(fp,Channels):
	A = {}
	f = open(fp)
	header = f.readline()
	for line in f:
		ls = line.strip().split('\t')
		A[ls[0]] = float(ls[4])
	f.close()
	return [A[c] for c in Channels]

def lookup_dict(fp):
	L = {}
	f = open(fp+'lookup.txt')
	for line in f:
		ls = line.strip().split()
		L[int(ls[0])] = int(ls[1])
	f.close()
	return L

strints = [str(x) for x in range(10)]
def parse_mcl_clusters(fp):
	C = [[]]
	f = open(fp)
	for line in f:
		if line[0] in strints:
			C[-1] += [int(x) for x in line.replace('$','').strip().split()[1:]]
		elif line[0] == ' ':
			C[-1] += [int(x) for x in line.replace('$','').strip().split()]
		if '$' in line:
			C.append([])
	f.close()
	return C

def cluster_stats(clusters,lookup,data,ids,thresh,min_size=100,c_id=None):
	C = []
	for c in clusters:
		if len(c) > min_size:
			idx = [ids[lookup[x]] for x in c]
			d = data[idx]
			idx_thresh = np.transpose([(d[:,i] > t) for i,t in enumerate(thresh)])
			alpha = idx_thresh.sum(0,dtype=np.float64)/d.shape[0]
			mu = np.zeros(len(alpha))
			sigma = np.zeros(len(alpha))
			for i in range(idx_thresh.shape[1]):
				if alpha[i] > 0:
					mu[i] = np.average(d[idx_thresh[:,i],i])
					sigma[i] = np.std(d[idx_thresh[:,i],i])
			if c_id:
				C.append((alpha,mu,sigma,set(idx),c_id))
			else:
				C.append((alpha,mu,sigma,set(idx)))
		else:
			break
	return C

def resolution_clustering(clusters,cluster_ids,sampled,kx=2):
	X = np.array([np.append(np.append(c[0],c[1]),c[2]) for c in clusters])
	n = X.shape[1]/3
	Xn = X/([np.average(X[:,:n])]*n + [np.average(X[:,n:2*n])]*n + [np.average(X[:,2*n:])]*n)
	C,e,nf = Pycluster.kcluster(Xn,len(clusters)/len(sampled)*kx)
	del Xn
	Cidx = defaultdict(list)
	for i,c in enumerate(C):
		Cidx[c].append(i)
	CStable = []
	for k,v in Cidx.items():
		members = set()
		for c in v:
			members.update(clusters[c][3])
		members = sorted(members)
		s = stability(members,cluster_ids)
		CStable.append((s,np.average(X[v],axis=0).reshape((3,X.shape[1]/3)),members))
	return CStable

def stability(ids,clusters):
	d = 0.
	C = [c[ids] for c in clusters]
	pool = Pool(8)
	results = pool.map(cluster_xor,((i,C) for i in xrange(1,len(ids))),chunksize=len(ids)//8)
	pool.close()
	pool.join()
	for r in results:
		d += r
	return d/(misc.comb(len(ids),2)*misc.comb(len(clusters),2))

def cluster_xor(args):
	i,C = args
	c0 = [c[:-i] for c in C]
	c1 = [c[i:] for c in C]
	cAND = [(c0[a] == c1[a]) for a in range(len(C))]
	return (distance.pdist(cAND,'hamming')*len(c0[0])).sum()

def stability_summary(fp,s1=.01,s2=.05):
	f = open(fp)
	IDs = []
	S = []
	Channels = None
	Alpha = []
	Mu = []
	Sigma = []
	Members = []
	for line in f:
		ls = line.strip().split('\t')
		if 'cluster:' in line:
			IDs.append(int(ls[0].split()[1]))
			S.append(float(ls[1].split()[1]))
		elif (line[0] == '\t') and (Channels == None):
			Channels = ls
		elif 'alpha' in line:
			Alpha.append([float(x) for x in ls[1:]])
		elif 'mu' in line:
			Mu.append([float(x) for x in ls[1:]])
		elif 'sigma' in line:
			Sigma.append([float(x) for x in ls[1:]])
		elif len(ls) > 99:
			Members.append([int(x) for x in ls])
	TotalCapture1 = set()
	TotalCapture2 = set()
	t = set()
	for i,m in enumerate(Members):
		if S[i] < s2:
			TotalCapture2.update(set(m))
			if S[i] < s1:
				TotalCapture1.update(set(m))
		t.update(set(m))
	return IDs,S,Channels,Alpha,Mu,Sigma,Members,len(TotalCapture1)/float(len(t)),len(TotalCapture2)/float(len(t))

def update_stable(R,C,I,A,M,S,l,s=.05):
	C.append(R[2])
	for i,id in enumerate(R[0]):
		if R[1][i] > s:
			break
		I.append('|'.join(l+[str(R[0][i]),str(R[1][i])]))
		A.append([str(x) for x in R[3][i]])
		M.append([str(x) for x in R[4][i]])
		S.append([str(x) for x in R[5][i]])
	return C,I,A,M,S

code_dir = '/ahg/regevdata/users/bcleary/CyTOF/code'
help_message = 'usage example: python StabilityAnalysis.py -i some_experiment/ -t alpha_thresh_file.tab'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hi:t:',["inputdir=","threshfile="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-h','--help'):
			print help_message
			sys.exit()
		elif opt in ('-i','--inputdir'):
			if arg[-1] == '/':
				inputdir = arg
			else:
				inputdir = arg+'/'
		elif opt in ('-t','--threshfile'):
			threshfile = arg
	Channels,Ids,Data = get_data(inputdir+'all_data.tab.filtered')
	Thresholds = get_alpha_thresholds(threshfile,Channels)
	FP = glob.glob(os.path.join(inputdir,'bootstrap','*/'))
	AllClusters = defaultdict(list)
	SampledIds = {}
	ClusterIds = defaultdict(list)
	for fp in FP:
		Lookup = lookup_dict(fp)
		ClusterFiles = glob.glob(os.path.join(fp,'out.*'))
		for cf in ClusterFiles:
			C = parse_mcl_clusters(cf)
			boot_id = cf[cf.index('bootstrap/')+10:cf.rfind('/')]
			resolution = cf[cf.rfind('.')+1:]
			C = cluster_stats(C,Lookup,Data,Ids,Thresholds,c_id=boot_id)
			Cid = np.random.random(len(Ids))
			for i,c in enumerate(C):
				AllClusters[resolution].append(c)
				Cid[list(c[3])] = i
			ClusterIds[resolution].append(Cid)
		SampledIds[boot_id] = set(Lookup.values())
	for k,v in AllClusters.iteritems():
		if len(v) > 1:
			S = sorted(resolution_clustering(v,ClusterIds[k],SampledIds),key=itemgetter(0))
			f = open(inputdir + 'bootstrap/cluster_stability.%s.tab' % k,'w')
			for i,s in enumerate(S):
				f.write('cluster: %d\tstability: %f\tnumber of members: %d\n\n' % (i,s[0],len(s[2])))
				f.write('\t'.join(['']+Channels) + '\n')
				for j,value in enumerate(['alpha','mu','sigma']):
					f.write('\t'.join([value] + [str(x) for x in s[1][j]]) + '\n')
				f.write('\n')
				f.write('\t'.join([str(x) for x in s[2]]) + '\n')
				f.write('\n\n')
			f.close()

