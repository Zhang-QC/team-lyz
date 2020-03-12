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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import zhangqc_sca as zs
import integrate_pymol as ip


def run1(pdb_id, n_seq, E_value):
	'''
	Run the program.
	
	Input:
		PDB_id: a string
		n_seq: integer
		E-value: integer

	Output:
		files saved in the folder /lyzsite/static
		the reference position in the aligned fasta file (integer)
	'''
	uniprot_id = util.get_uniprot_id(pdb_id)
	my_protein = util.Protein(pdb_id, uniprot_id)
	similar = util.get_similar(my_protein.name, n_seq)
	print("Completed finding similar sequences, generationg Multiple Sequence Alignment:")
	#record = SeqRecord(Seq(my_protein.sequence, IUPAC.protein), 
	#	id = my_protein.uniprot_id, name=my_protein.name)
	#with open('./lyzsite/static/reference.fasta', 'w') as output_handle:
	#	SeqIO.write(record, output_handle, 'fasta')
	if my_protein.uniprot_id not in similar:
		similar.append(my_protein.uniprot_id)
	my_msa = util.create_MSA(similar)
	records = list(SeqIO.parse("./lyzsite/static/aligned.fasta", "fasta"))
	for index, seq in enumerate(records):
		if seq.id == uniprot_id:
			return index


def run2(pdb_id, ref_id, fasta_filename):
	print('Performing the Statistical Coupling Analysis')
	output_name = 'sca_result.db'
	zs.perform_calculations(fasta_filename,\
	 pdb_id, ref_id, output_name)


def run3():
	print('Generating output plots')
	if not os.path.exists('Outputs/'): 
		os.makedirs('Outputs/')
	db, Dseq, Dsca, Dsect, listS, ind = zs.process_output('Outputs/aligned.db')
	zs.image_pairwise(Dseq, Dsca, Dsect, listS, ind)
	zs.image_conservation(Dseq, Dsca, Dsect, listS, ind)
	zs.image_matrix(Dseq, Dsca, Dsect, listS, ind)
	#zs.image_eigenvalues(Dseq, Dsca, Dsect, listS, ind)
	os.chdir('..')


def run4(pdb_id):
	print('Generating pyMol graph')
	#os.chdir('./lyzsite/static')
	ip.create('model.py','5', pdb_id)
	os.system('pymol model.py')


def run_all(pdb_id, n_seq, E_value):
	index = run1(pdb_id, n_seq, E_value)
	run2(pdb_id, index, "./lyzsite/static/aligned.fasta")
	run3()
	run4(pdb_id)

