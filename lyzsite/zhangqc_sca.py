from __future__ import division

import os
import time
import matplotlib.pyplot as plt
import numpy as np
import copy
import scipy.cluster.hierarchy as sch
from scipy.stats import scoreatpercentile
import matplotlib.image as mpimg
from IPython.display import display
from IPython.display import Image
from Bio.Seq import Seq
from Bio import motifs
import colorsys
import pySCA.scaTools as sca
import mpld3
import pickle
from optparse import OptionParser


def perform_calculations(fasta_name, pdb_id, ref_id, output_name):
	#os.chdir('./pySCA')
	cwd = os.getcwd()
	print('***', cwd)
	l1 = 'python3 ./pySCA/scaProcessMSA.py ' + fasta_name + ' -i ' + str(ref_id)\
	 + ' -t' #+ '--output ' + output_name 
	l2 = 'python3 ./pySCA/scaCore.py Outputs/aligned.db'
	l3 = 'python3 ./pySCA/scaSectorID.py Outputs/aligned.db'
	print(l1)
	os.system(l1) 
	os.system(l2)
	#os.system(l3)
	#os.chdir('..')


def process_output(output_name):
	cwd = os.getcwd()
	print('$$$', cwd)
	with open(output_name, 'rb') as handle: 
		db = pickle.load(handle) 
	Dseq = db['sequence']
	Dsca = db['sca']
	Dsect = None #db['sector']
	listS = [Dsca['simMat'][i,j] for i in range(Dsca['simMat'].shape[0]) \
			 for j in range(i+1, Dsca['simMat'].shape[1])]
	Z = sch.linkage(Dsca['simMat'],method = 'complete', metric = 'cityblock')
	R = sch.dendrogram(Z, no_plot = True)
	ind = map(int, R['ivl'])
	return db, Dseq, Dsca, Dsect, listS, ind


def image_pairwise(Dseq, Dsca, Dsect, listS, ind):
	plt.rcParams['figure.figsize'] = 9, 4
	plt.subplot(121)
	plt.hist(listS, Dseq['Npos'])
	plt.xlabel('Pairwise sequence identities', fontsize=14)
	plt.ylabel('Number', fontsize=14)
	plt.subplot(122)
	plt.imshow(Dsca['simMat'], vmin=0, vmax=1) 
	plt.colorbar()
	plt.savefig('lyzsite/static/1.png')


def image_conservation(Dseq, Dsca, Dsect, listS, ind):
	fig, axs = plt.subplots(1,1, figsize=(9,4))
	xvals = [i+1 for i in range(len(Dsca['Di']))]
	xticks = [0,45,95,144]
	plt.bar(xvals,Dsca['Di'], color='k')
	plt.tick_params(labelsize=11); plt.grid()
	axs.set_xticks(xticks);
	labels = [Dseq['ats'][k] for k in xticks]
	axs.set_xticklabels(labels);
	plt.xlabel('Amino acid position', fontsize=18); plt.ylabel('Di', fontsize=18);
	plt.savefig('lyzsite/static/2.png')


def image_matrix(Dseq, Dsca, Dsect, listS, ind):
	plt.rcParams['figure.figsize'] = 13, 8
	plt.imshow(Dsca['Csca'], vmin=0, vmax=1.4,interpolation='none',\
			   aspect='equal')
	plt.savefig('lyzsite/static/3.png')


'''
def image_eigenvalues(Dseq, Dsca, Dsect, listS, ind):
	plt.rcParams['figure.figsize'] = 9, 4
	hist0, bins = np.histogram(Dsca['Lrand'].flatten(), bins=Dseq['Npos'], \
						   range=(0,Dsect['Lsca'].max()))
	hist1, bins = np.histogram(Dsect['Lsca'], bins=Dseq['Npos'], \
							   range=(0,Dsect['Lsca'].max()))
	plt.bar(bins[:-1], hist1, np.diff(bins),color='k')
	plt.plot(bins[:-1], hist0/Dsca['Ntrials'], 'r', linewidth=3)
	plt.tick_params(labelsize=11)
	plt.xlabel('Eigenvalues', fontsize=18); plt.ylabel('Numbers', fontsize=18)
	plt.savefig('eigenvalues.jpg')
'''

def find_sectors(Dseq, Dsca, Dsect, listS, ind):
	pass





