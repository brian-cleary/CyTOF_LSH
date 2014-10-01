import numpy as np

def cluster_zscores(C,X,L):
	all_idx = []
	for c in C:
		all_idx += c
	Xavg = np.average(X[all_idx],axis=0)
	Xstd = np.std(X[all_idx],axis=0)
	Z = []
	for c in C:
		Cavg = np.average(X[c],axis=0)
		z = (Cavg - Xavg)/Xstd
		Z.append([(L[i],z[i]) for i in np.argsort(np.abs(z))[::-1]])
	return Z

def enrichment_by_experiment(Clusters,Xlabels,Experiments,n=100):
	#Experiments looks like: {'LN.131030': 0,'somethingelse': 1,...}
	Counts = np.zeros((len(Clusters),len(Experiments)))
	newClusters = []
	for i,c in enumerate(Clusters):
		nc = []
		for j,l in enumerate(Xlabels[c]):
			if l in Experiments:
				Counts[i,Experiments[l]] += 1
				nc.append(c[j])
		newClusters.append(nc)
	Fracs = Counts/Counts.sum(0)
	i = np.argsort(Fracs.sum(1))[:-n-1:-1]
	i = [j for j in i if Counts[j].sum() > 0]
	return Counts[i],i,[newClusters[j] for j in i]

Exp = {}
for l in set(Xl):
	if 'CNS' in l:
		Exp[l] = len(Exp)

Alpha = []
Mu = []
Sigma = []
Source = []
for c in C[:2000]:
	a = []
	m = []
	s = []
	sc = []
	for j in range(X.shape[1]):
		x = X[np.nonzero(X[c,j])[0],j]
		a.append(len(x)/float(len(c)))
		m.append(x.sum()/len(c))
		if len(x) > 0:
			s.append(np.std(x))
		else:
			s.append(0)
	nx = [N[x] for x in c]
	Source.append([nx.count(n)/x for n,x in Ns])
	Alpha.append(a)
	Mu.append(m)
	Sigma.append(s)