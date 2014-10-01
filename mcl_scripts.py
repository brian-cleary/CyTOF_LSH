#!/usr/bin/env python

import sys,getopt
import numpy as np
from scipy import cluster
from scipy.spatial import distance
from collections import defaultdict

def cluster_histograms(cl):
	DateHist = defaultdict(int)
	OtherHist = defaultdict(int)
	for l in cl:
		DateHist[l[:l.index('.')]] += 1
		OtherHist[l[l.index('.')+1:l.index('|')]] += 1
	return DateHist,OtherHist

def plot_hists(c,L,Xd,N,Ns):
	ticks = [[int(t*(5-i)/5) for i in range(5)] for t in Xd.max(0)]
	for i in range(len(L)):
		pylab.subplot(6,6,i+1)
		pylab.hist(Xd[c,i],log=True)
		pylab.xticks(ticks[i])
		pylab.title(L[i])
	pylab.subplot(6,6,i+2)
	nx = [N[x] for x in c]
	pylab.plot([nx.count(n)/x for n,x in Ns],'*')
	pylab.title('source')
	pylab.xticks(np.arange(len(Ns)),[ns[0] for ns in Ns],rotation='vertical',size=8)
	pylab.subplots_adjust(hspace=.8,wspace=.22,left=0.06,bottom=0.10,right=0.95)
	pylab.show()

def dendro_clusters(C,X,ct=.2):
	CX = [np.average(X[c],axis=0) for c in C]
	D = distance.pdist(CX,'minkowski',5)
	L = hierarchy.linkage(D,method='complete',metric='minkowski')
	hierarchy.dendrogram(L,color_threshold=ct,labels=np.arange(len(C)),leaf_font_size=9)
	pylab.show()

def return_clusters(C,X,D=None,ct=0.42,clust_meth='complete'):
	CX = [np.average(X[c],axis=0) for c in C]
	D = distance.pdist(CX,'minkowski',5)
	L = hierarchy.linkage(D,method=clust_meth,metric='minkowski')
	CC = hierarchy.fcluster(L,ct,criterion='distance')
	CD = defaultdict(list)
	for i,x in enumerate(CC):
		CD[x].append(i)
	CC = []
	for v in CD.values():
		x = []
		for c in v:
			x += C[c]
		CC.append((len(x),x))
	CC.sort(reverse=True)
	return [c[1] for c in CC]

import numpy as np
from collections import defaultdict
from scipy.spatial import distance
from scipy.cluster import hierarchy
def cluster_divergence(CX,clust_meth='complete'):
	#CX = [np.average(X[c],axis=0) for c in C]
	print 'doing distance'
	D = distance.squareform(distance.pdist(CX,'minkowski',5))
	print 'doing linkage'
	L = hierarchy.linkage(D,method=clust_meth,metric='minkowski')
	divergence = []
	for ct in np.linspace(.25,1.25):
		CC = hierarchy.fcluster(L,ct,criterion='distance')
		CD = defaultdict(list)
		for i,x in enumerate(CC):
			CD[x].append(i)
		CD = CD.values()
		d = []
		for c in CD:
			if len(c) > 1:
				x = []
				for i in range(len(c)):
					for j in range(i+1,len(c)):
						x.append(D[i,j])
				d.append(max(x))
		divergence.append((ct,np.average(d)))
		print divergence[-1]
	np.save('cluster_diameters.npy',divergence)
			

f = open('dat/all_data.tab')
L = f.readline().strip().split('\t')[1:]
X = []
N = []
Ns = defaultdict(float)
for line in f:
	ls = line.strip().split('\t')
	n = ls[0].split('|')[0]
	N.append(n)
	Ns[n] += 1
	X.append([float(x) for x in ls[1:]])

f.close()
X = np.array(X)
Ns = sorted(Ns.items())

FP = glob.glob(os.path.join('*','all_data.tab.norm'))
L = None
for fp in FP:
	f = open(fp)
	l = f.readline().strip().split('\t')
	if L:
		L = (L & set(l))
	else:
		L = set(l)

L = list(L)
Lf = []
for l in L:
	add = True
	for fil in ['DNA','Beads','bead','Xe131Di','Event','Time','Ba138Di','Cs133Di','BCKG190Di','Viability']:
		if fil in l:
			add = False
			break
	if add:
		Lf.append(l)

del Lf[Lf.index('ID')]

g = open('all_clusters.I16.tab','w')
g.write('\t'.join(['ID']+Lf)+'\n')
ga = open('all_clusters.I16.alpha.tab','w')
ga.write('\t'.join(['ID']+Lf)+'\n')
gm = open('all_clusters.I16.mu.tab','w')
gm.write('\t'.join(['ID']+Lf)+'\n')
for expr in ['Treg','Th0','Th1','Th2','B623','T16','T36']:
	f = open('131218_72h_%s/all_data.tab'%expr)
	lf = f.readline().strip().split('\t')[1:]
	Li = [lf.index(l) for l in Lf]
	Y = [[float(x) for x in line.strip().split('\t')[1:]] for line in f]
	Y = np.array(Y)
	Y = Y[:,Li]
	C = parse_mcl_clusters('131218_72h_%s/out.all_data.hash_trimmed.mci.I16'%expr)
	C = [c for c in C if len(c) > 9]
	Yavg = [np.average(Y[c],axis=0) for c in C]
	for i,y in enumerate(Yavg):
		g.write('\t'.join(['%s.%d' % (expr,i)] + [str(x) for x in y])+'\n')
	Yalpha = [(Y[c] > 0).sum(0)/float(len(c)) for c in C]
	for i,y in enumerate(Yalpha):
		ga.write('\t'.join(['%s.%d' % (expr,i)] + [str(x) for x in y])+'\n')
	Ymu = np.array([Y[c].sum(0)/(Y[c] > 0).sum(0) for c in C])
	Ymu[np.isnan(Ymu)] = 0
	for i,y in enumerate(Ymu):
		gm.write('\t'.join(['%s.%d' % (expr,i)] + [str(x) for x in y])+'\n')
	Ysigma = []
	for c in C:
		y = Y[c]
		sig = []
		for yt in y.T:
			sig.append(np.std(yt[np.nonzero(yt)]))
		Ysigma.append(sig)
	Ysigma = np.array(Ysigma)
	Ysigma[np.isnan(Ysigma)] = 0

f = open('all_clusters.norm.I16.alpha.tab')
g = open('all_clusters.norm.I16.alpha.zscores.tab','w')
g.write(f.readline())
X = []
Xl = []
for line in f:
	ls = line.strip().split('\t')
	Xl.append(ls[0])
	X.append([float(x) for x in ls[1:]])

X = np.array(X)
X = (X - np.average(X,0))/np.std(X,0)
X[np.isnan(X)] = 0
for i,x in enumerate(X):
	g.write('\t'.join([Xl[i]]+[str(y) for y in x])+'\n')

# mcl commands:
# create level 0 network from positive correlations
#	$ mcxarray -data wavelet.0.tab -co 0.6 -skipc 1 -o wavelet.0.60.mci -write-tab wavelet.0.dict
# find appropriate correlation level
#	$ mcx query -imx wavelet.0.60.mci --vary-correlation
# alter the network to remove correlations below 0.9
#	$ mcx alter -imx wavelet.0.60.mci -tf 'gq(0.9), add(-0.9)' -o wavelet.0.90.mci
#
# create level 3 network from anticorrelations
#	$ mcxarray -data wavelet.3.tab -co 0.2 -skipc 1 -tf 'mul(-1)' -o wavelet.3.20.mci
# find appropriate correlation level
#	$ mcx query -imx wavelet.3.20.mci --vary-correlation
# alter the network to remove correlations below 0.52
#	$ mcx alter -imx wavelet.3.20.mci -tf 'gq(0.52)' -o wavelet.3.52.mci
#
# find the interacting pairs from level 3
#	>>> X0 = get_interactions(path+'wavelet.3.52.mci')
# merge level 3 interactions with level 0 network
#	>>> merge_interactions(path+'wavelet.0.90.mci',path+'wavelet.merged.mci',X0)
#
# cluster the merged network
#	$ mcl wavelet.merged.mci -I 1.2
#	$ mcl wavelet.merged.mci -I 1.4
# write labels on cluster nodes
#	$ mcxdump -icl out.wavelet.merged.mci.I12 -o dump.wavelet.merged.mci.I12 -tabr wavelet.0.dict
#	$ mcxdump -icl out.wavelet.merged.mci.I14 -o dump.wavelet.merged.mci.I14 -tabr wavelet.0.dict

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

def entropy(v,d):
	x = np.zeros(len(d))
	for i in v:
		x[d[i]] += 1
	p = x/x.sum()
	p = p[np.nonzero(p)]
	return -sum(p*np.log(p))

g1 = open('cluster_sizes.tab','w')
g1.write('\t'.join(['Experiment','Inflation','Number of clusters','25th pct','Mean','Median','75th pct','Max'])+'\n')
g2 = open('cluster_accuracy.tab','w')
g2.write('\t'.join(['Experiment','Inflation','ClusterId','Size','Diversity'])+'\n')
E = ['B623','T16','T36','Th0','Th1','Th2','Treg']
I = ['I13','I16','I25']
for e in E:
	f = open(e+'/all_data.tab')
	line = f.readline()
	D = {}
	for line in f:
		ls = line.strip().split('\t')
		D[int(ls[0])] = float(ls[1])
	for i in I:
		C = parse_mcl_clusters('%s/out.all_data.hash_trimmed.mci.%s' % (e,i))
		C = [c for c in C if len(c) > 9]
		cs = [len(c) for c in C]
		g1.write('\t'.join([e,i] + [str(s) for s in [len(C),np.percentile(cs,25),np.average(cs),np.median(cs),np.percentile(cs,75),max(cs)]])+'\n')
		for j,c in enumerate(C):
			ent = entropy(c,D)
			g2.write('\t'.join([e,i] + [str(s) for s in [j,len(c),np.exp(ent)]])+'\n')

def parse_cytof_txt(fp):
	f = open(fp)
	L = f.read()
	f.close()
	L = [l.split('\t') for l in L.split('\r')]
	X = []
	for l in L[1:]:
		X.append([int(x) for x in l if x])
	return L[0],np.array(X)

def parse_cytof_txt(fp,spt='\t'):
	f = open(fp)
	line = f.readline()
	L = f.readline().strip().split(spt)
	X = [[float(x) for x in line.strip().split(spt)] for line in f]
	f.close()
	return L,np.array(X)

def gateX(X,L,dnaG,viadnaG):
	m = (viadnaG[1][1] - viadnaG[0][1])/(viadnaG[1][0] - viadnaG[0][0])
	line_eqn = [-m*viadnaG[0][0]-viadnaG[0][1],m,-1]
	via_idx = L.index('Viability')
	dna_idx = [i for i,l in enumerate(L) if 'DNA' in l]
	Y = np.where((X[:,dna_idx[0]] > dnaG[0]) & (X[:,dna_idx[1]] > dnaG[1]) & (np.dot(X[:,[dna_idx[0],via_idx]],line_eqn[1:]) > line_eqn[0]) & (X[:,via_idx] < viadnaG[1][1]))
	return X[Y[0]]

def plot_controls(X,L):
	time_idx = L.index('Time')
	via_idx = L.index('Viability')
	dna_idx = [i for i,l in enumerate(L) if 'DNA' in l]
	pylab.subplot(311)
	pylab.loglog(X[:,dna_idx[0]],X[:,dna_idx[1]],'.')
	pylab.subplot(312)
	pylab.loglog(X[:,dna_idx[0]],X[:,via_idx],'.')
	pylab.subplot(313)
	pylab.loglog(X[:,dna_idx[0]],X[:,time_idx],'.')
	pylab.show()

def set_of_labels(R):
	Ls = set(R[0][0])
	for r in R[1:]:
		Ls = set(Ls & set(r[0]))
	return Ls

def merged_matrix(R,row_labels,col_labels):
	X = []
	for i,r in enumerate(R):
		Lj = [r[0].index(l) for l in col_labels]
		for j,x in enumerate(r[1]):
			xj = [x[lj] for lj in Lj]
			if sum(xj) > 0:
				X.append([row_labels[i]+'|'+str(j)] + xj)
	return X

def merge_and_normalize(R,row_labels,col_labels):
	X = []
	for i,r in enumerate(R):
		j,n = normalize_matrix(r[0],r[1],col_labels)
		for k in range(len(j)):
			X.append([row_labels[i]+'|'+str(j[k])] + list(n[k]))
	return X

def normalize_matrix(L,X,col_labels,p=5):
	Lj = [L.index(l) for l in col_labels]
	y = []
	j = []
	for i,x in enumerate(X):
		if x:
			y.append(x)
			j.append(i)
	y = np.array(y)[:,Lj]
	y = np.arcsinh(y)
	return j,y/(np.sum(y**p,axis=-1)**(1./p))[:,np.newaxis]

def read_block(f,d={}):
	for _ in range(500):
		line = f.readline()
		ls = line.strip().split('\t')
		d[int(ls[0])] = ls[1:]
	return d

def merge_distances(op,max_edges,min_corr):
	FP = glob.glob(os.path.join(op,'*.distances.tab'))
	F = [open(fp) for fp in FP]
	D = [read_block(f) for f in F]
	
