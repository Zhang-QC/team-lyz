'''
***********
Jason Zhang
***********
'''

import bs4
import requests
import util
import tuo
#import jenny


def get_uniprot_id(pdb_id):
	'''
	Input a Protein Data Base ID code, produce the associated UniProt ID code.
	Return None if the PDB ID is not valid.

	Input:
		pdb_id: a string

	Output:
		A string
	'''
	pdb_url = "http://www.rcsb.org/structure/" + pdb_id
	r = requests.get(pdb_url)
	text = r.text.encode('iso-8859-1')
	soup = bs4.BeautifulSoup(text, "html5lib")

	url_tag = soup.find_all('a')
	for potential_url in url_tag:
		if potential_url.has_attr('href'):
			check_url = potential_url['href']
			if check_url[0] != "#":
				if check_url[:31] == "http://www.uniprot.org/uniprot/":
					return check_url[31:]
	print("Your PDB ID is not valid.")
	return None


def is_sequence(seq):
	'''
	Check if a sequence is a valid amino acid sequence. A valid AA sequence 
	must have all uppercase letters that corresponds to the 20 one-letter amino
	acid codes.
	'''
	for i in seq:
		if i not in util.AMINO_ACIDS:
			return False
	return True


FASTA_EXAMPLE = '''
>sp|P07830|ACT1_DICDI Major actin OS=Dictyostelium discoideum OX=44689 GN=act1 PE=1 SV=2
MDGEDVQALVIDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHTGVMVGMGQKDSYVGDEAQ
SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTEAPLNPKANREKM
TQIMFETFNTPAMYVAIQAVLSLYASGRTTGIVMDSGDGVSHTVPIYEGYALPHAILRLD
LAGRDLTDYMMKILTERGYSFTTTAEREIVRDIKEKLAYVALDFEAEMQTAASSSALEKS
YELPDGQVITIGNERFRCPEALFQPSFLGMESAGIHETTYNSIMKCDVDIRKDLYGNVVL
SGGTTMFPGIADRMNKELTALAPSTMKIKIIAPPERKYSVWIGGSILASLSTFQQMWISK
EEYDESGPSIVHRKCF
'''

def parse_fasta_header(st):
	'''
	Taking in a FASTA header string (always starting with '>'), return
	the information contained in tuple form. Some of the values, may
	be missing and will be returned as None.

	The documentation for FASTA file is taken from:
	https://www.uniprot.org/help/fasta-headers
	
	Input:
		st: a string

	Output:
		a tuple containing:
			db: a string (either 'sp' for UniProtKB/Swiss-Prot
			 or 'tr' for UniProtKB/TrEMBL)
			identifier: a string containing the UniProt ID
			entry_name: the name of the UniProt entry
			protein_name: the name of the protein (annotated)
			organism_name (OS): the scientific name of the organism
			organism_identifer (OX): the unique identifier of the organism
			gene_name (GN): the gene name of the UniProtKB entry
			protein_existence (PE): the value describing the existence
			 of the protein
			sequence_version (SV): the version number of the sequence
	'''	
	l = st.split('|')
	db = l[0][1:]
	identifier = l[1]
	l = l[2].split(' ')
	dic = {'OS': [], 'OX': [], 'GN': [], 'PE': [], 'SV': []}
	entry_name = l[0]
	protein_name = []
	i = 1
	while '=' not in l[i]:
		protein_name.append(l[i])
		i += 1
	protein_name = ' '.join(protein_name)
	l = l[i:]
	for j in l:
		if '=' in j:
			s = j[0: 2]
			dic[s].append(j[3:])
		else:
			dic[s].append(j)
	for k in dic:
		if dic[k] == []:
			dic[k] = None
		else:
			dic[k] = ' '.join(dic[k])
	return db, identifier, entry_name, protein_name, dic


def parse_fasta(st):
	'''
	Taking in a FASTA string, return the name species and sequence of the 
	associated protein.
	
	Input:
		st: a string from a FASTA file

	Output:
		header: string
		name: string
		species: string
		sequence: string
	'''
	l = st.split('\n')
	while '' in l:
		l.remove('')
	seq = ''	
	for i in l:
		if i[0] == '>':
			header = i
		else:
			seq += i
	db, identifier, entry_name, protein_name, dic = parse_fasta_header(header)
	return header, protein_name, dic['OS'], seq
	