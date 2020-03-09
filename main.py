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
import jenny
import os
from Bio import SeqIO
from Bio import AlignIO
import zhangqc_sca as zs
import integrate_pymol as ip


def run(pdb_id, n_seq, E_value):
	'''
	Run the program.
	
	Input:
		PDB_id: a string
		n_seq: integer
		E-value: integer

	Output:
		files saved in the folder /lyzsite/static
	'''
	uniprot_id = util.get_uniprot_id(pdb_id)
	my_protein = util.Protein(pdb_id, uniprot_id)
	similar = util.get_similar(my_protein.name, n_seq)
	my_msa = util.create_MSA(similar)


def run2(pdb_id, my_protein, fasta_filename):
	#os.chdir('pySCA')
	output_name = './lyzsite/static/output.db'
	zs.perform_calculations(fasta_filename,\
	 pdb_id, my_protein.species, output_name)


def run3():
	os.chdir('..')
	zs.process_output(output_name)
	zs.image_pairwise(Dseq, Dsca, Dsect, listS, ind)
	zs.image_conservation(Dseq, Dsca, Dsect, listS, ind)
	zs.image_matrix(Dseq, Dsca, Dsect, listS, ind)
	zs.image_eigenvalues(Dseq, Dsca, Dsect, listS, ind)


def run4(pdb_id):
	#os.chdir('./lyzsite/static')
	ip.create('model.py','5', pdb_id)
	os.system('pymol model.py')

