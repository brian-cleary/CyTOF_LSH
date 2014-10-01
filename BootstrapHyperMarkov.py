#!/usr/bin/env python

import sys, getopt, os
import random
import numpy as np

def get_all_data(fp):
	f = open(fp)
	header = f.readline()
	Lines = f.readlines()
	return header,Lines

def sample_points(header,Lines,k,fp):
	f = open(fp,'w')
	f1 = open(fp[:fp.rfind('/')]+'/lookup.txt','w')
	f.write(header)
	for i,line in enumerate(random.sample(Lines,k)):
		f.write(line)
		f1.write('%d %s\n' % (i,line.split()[0]))
	f.close()
	f1.close()

def tune_hash_size(num_rows,num_cols,d=5,sim=0.8,p=0.95):
	# assume the tree will walk np.log2(num_rows/num_cols/d) steps
	s = np.log2(num_rows/num_cols/d)
	h = np.log(1-p)/np.log(1 - (1 - np.arccos(sim)/np.pi )**s)
	return int(np.round(h))

code_dir = '/ahg/regevdata/users/bcleary/CyTOF/code'
help_message = 'usage example: python BootstrapHyperMarkov.py -i some_experiment/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'hi:r:s:',["inputdir=","rounds=","sample="])
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
		elif opt in ('-r','--rounds'):
			num_bootstraps = int(arg)
		elif opt in ('-s','--sample'):
			sample_fraction = float(arg)
	# assuming parsed, normalized, etc. data...
	H,L = get_all_data(inputdir + 'all_data.tab.norm.filtered')
	sample_fraction = int(np.round(len(L)/sample_fraction))
	inputdir += 'bootstrap/'
	hash_size = tune_hash_size(sample_fraction,len(H.split('\t'))-1)
	for b in range(num_bootstraps):
		current_dir = '%s/%d/' % (inputdir,b)
		os.system('mkdir %s' % current_dir)
		os.system('mkdir %s/Logs' % current_dir)
		os.system('cp %s/*.q %s' % (code_dir,current_dir))
		os.system('sed -i "s/\/project\/home\//%s/g" %s/*.q' % (current_dir.replace('/','\/'),current_dir))
		os.system('sed -i "s/HashDist\[1-50\]/HashDist\[1-%d\]/g" %s/HashDist_ArrayJob.q' % (hash_size,current_dir))
		os.system('sed -i "s/JoinDist\[1-150\]/JoinDist\[1-%d\]/g" %s/JoinDist_ArrayJob.q' % (int(np.ceil(sample_fraction/1000.)),current_dir))
		sample_points(H,L,sample_fraction,current_dir+'sampled_data.tab.norm.filtered')
		os.system('mkdir %s/distances' % current_dir)
		os.system('sed -i "s/HashDist/HashDist%d/g" %s/HashDist_ArrayJob.q' % (b,current_dir))
		os.system('bsub < %s/HashDist_ArrayJob.q' % current_dir)
		os.system('mkdir %s/joins' % current_dir)
		os.system('sed -i "s/JoinDist/JoinDist%d/g" %s/JoinDist_ArrayJob.q' % (b,current_dir))
		os.system('bsub -w "ended(HashDist%d)" < %s/JoinDist_ArrayJob.q' % (b,current_dir))
		os.system('bsub -w "ended(JoinDist%d)" -q hour -o %s/Logs/rm.out -e %s/Logs/rm.err rm -r %s/distances' % (b,current_dir,current_dir,current_dir))
		os.system('sed -i "s/MCLCluster/MCLCluster%d/g" %s/MCLCluster_Job.q' % (b,current_dir))
		os.system('bsub -w "ended(JoinDist%d)" < %s/MCLCluster_Job.q' % (b,current_dir))
		os.system('bsub -w "ended(MCLCluster%d)" -q hour -o %s/Logs/rm.out -e %s/Logs/rm.err rm -r %s/joins %s/sampled_data.tab.norm.filtered %s/sampled_data.hash_trimmed.mci' % (b,current_dir,current_dir,current_dir,current_dir,current_dir))
		# ids need to be reconciled
		
	# estimate cluster stability within given resolution:
	#	flat clustering using complete distance from KL or alpha/sigma/mu vectors
	#	take only those clusters with >n elements, merge elements into "stable" cluster
	#for r in resolutions: (highest res first)
	#	if (cluster is big enough) and (stable):	??care about only writing clusters we haven't seen??
	#		write it down