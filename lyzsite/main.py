import bs4
import requests
import re
import urllib3
import certifi
import pa2util
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment
import util
import zhangqc
import os
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import zhangqc_sca as zs
import integrate_pymol as ip


def run1(pdb_id, n_seq, l_value):
	'''
	Run the program.
	
	Input:
		PDB_id: a string
		n_seq: a string
		l_value: a string

	Output:
		files saved in the folder /lyzsite/static
		the reference position in the aligned fasta file (integer)
	'''
	uniprot_id = util.get_uniprot_id(pdb_id)
	my_protein = util.Protein(pdb_id, uniprot_id)
	similar = util.get_similar(my_protein.name, n_seq)
	n_seq = int(n_seq)
	l_value = float(l_value)
	print("Completed finding similar sequences, generationg Multiple Sequence\
	 Alignment:")
	if my_protein.uniprot_id not in similar:
		similar.append(my_protein.uniprot_id)
	max_len = my_protein.length * l_value
	my_msa = util.create_MSA(similar, max_len)
	records = list(SeqIO.parse("./static/aligned.fasta", "fasta"))
	for index, seq in enumerate(records):
		if seq.id == uniprot_id:
			return index


def run2(pdb_id, ref_id, fasta_filename):
	'''
	Run the program to calculate statistical coupling analysis
	
	Input:
		PDB_id: a string
		ref_id: a string
		fasta_filename: the name of the MSA fasta files

	Output: None
	'''
	print('Performing the Statistical Coupling Analysis')
	output_name = 'sca_result.db'
	zs.perform_calculations(fasta_filename,\
	 pdb_id, ref_id, output_name)


def run3():
	'''
	Run the program to generate plots for output
	
	Input: 
		None

	Output: 
		None
	'''
	print('Generating output plots')
	if not os.path.exists('Outputs/'): 
		os.makedirs('Outputs/')
	cwd = os.getcwd()
	print('&&&', cwd)
	db, Dseq, Dsca, Dsect, listS, ind = zs.process_output('Outputs/aligned.db')
	zs.image_pairwise(Dseq, Dsca, Dsect, listS, ind)
	zs.image_conservation(Dseq, Dsca, Dsect, listS, ind)
	zs.image_matrix(Dseq, Dsca, Dsect, listS, ind)
	return db, Dseq, Dsca, Dsect, listS, ind


def run4(pdb_id, col_lst):
	'''
	Run the program to generate pyMol graphs
	
	Input: 
		None

	Output: 
		None
	'''
	print('Generating pyMol graph')
	os.chdir('static')
	ip.create('model.py','4','5', pdb_id, col_lst)
	os.system('pymol model.py')
	os.chdir('..')
	

def run_all(pdb_id, n_seq, E_value):
	'''
	Run the program to perform all the functions above using
	the users' inputs
	
	Input:
		PDB_id: a string
		n_seq: integer
		E-value: integer

	Output: 
		None
	'''
	index = run1(pdb_id, n_seq, E_value)
	run2(pdb_id, index, "./static/aligned.fasta")
	db, Dseq, Dsca, Dsect, listS, ind = run3()
	run4(pdb_id, Dsca['Di'])
