import bs4
import requests
import re
import urllib3
import certifi
import pa2util
import os
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align.Applications import ClustalOmegaCommandline as coc


AMINO_ACIDS = {"A": "alanine", "R": "arginine", "N": "asparagine", 
"D": "aspartate", "C": "cysteine", "E": "glutamate", "Q": "glutamine",
"G": "glycine", "H": "histidine", "I": "isoleucine", "L": "leucine",
"K": "lysine", "M": "methionine", "F": "phenylalanine", "P": "proline",
"S": "serine", "T": "threonine", "W": "tryptophan", "Y": "tyrosine",
"V": "valine"}


FASTA_EXAMPLE = '''
>sp|P07830|ACT1_DICDI Major actin OS=Dictyostelium discoideum OX=44689 \
GN=act1 PE=1 SV=2
MDGEDVQALVIDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHTGVMVGMGQKDSYVGDEAQ
SKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTEAPLNPKANREKM
TQIMFETFNTPAMYVAIQAVLSLYASGRTTGIVMDSGDGVSHTVPIYEGYALPHAILRLD
LAGRDLTDYMMKILTERGYSFTTTAEREIVRDIKEKLAYVALDFEAEMQTAASSSALEKS
YELPDGQVITIGNERFRCPEALFQPSFLGMESAGIHETTYNSIMKCDVDIRKDLYGNVVL
SGGTTMFPGIADRMNKELTALAPSTMKIKIIAPPERKYSVWIGGSILASLSTFQQMWISK
EEYDESGPSIVHRKCF
'''


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
	soup = bs4.BeautifulSoup(text, "html.parser")
	url_tag = soup.find_all('a')
	for potential_url in url_tag:
		if potential_url.has_attr('href'):
			check_url = potential_url['href']
			if check_url[0] != "#":
				if check_url[:31] == "http://www.uniprot.org/uniprot/":
					return check_url[31:]
	print("Your PDB ID is not valid.")
	return None


def get_fasta(uniprot_id):
	'''
	Get Fasta sequence with uniprot id of a protein

	Input:
		uniprot_id: a string

	Output:
		text: fasta sequence string
	'''
	uniprot_url = "https://www.uniprot.org/uniprot/" + uniprot_id + ".fasta"
	r = requests.get(uniprot_url)
	text = str(r.text.encode())[2:].replace('\\n', '\n')
	return text


def read_fasta(fast_str):
	'''
	Process a fasta sequence string by getting rid of undesired characters

	Input:
		fasta_str: fasta sequence

	Output: 
		total_list: a list of subsequences of a fasta string
	'''
	strin = fast_str.split('\n')
	substr = strin[0]
	sub_list = re.findall('([A-Z][0-9]{5}),([A-Z0-9]{4}_[A-Z]{5}),(OS=[A-Z]\
		[a-z\s]+),(OX=[0-9]{4}),(GN=[A-Za-z0-9]{5}),(PE=[0-9]),(SV=[0-9])',\
		 substr)
	total_list = sub_list
	for sub_list in strin[1:]:
		if sub_list != ' ':
			total_list.append(sub_list)
	return total_list
	

def code_search(url):
	'''
	Collect all the uniprot ids of proteins that are existing on the webpage

	Input:
		url: a given url

	Output: 
		lst_codes: a list of uniprot id
	'''
	r = requests.get(url)
	text = r.text.encode()
	soup = bs4.BeautifulSoup(text,"html.parser")
	table = soup.find("table")
	trs = table.find("tbody").find_all('tr')
	lst_codes = []
	for row in trs:
		lst_codes.append(row.attrs["id"])
	return lst_codes


def is_sequence(seq):
	'''
	Check if a sequence is a valid amino acid sequence. A valid AA sequence 
	must have all uppercase letters that corresponds to the 20 one-letter 
	amino acid codes.
	'''
	for i in seq:
		if i not in util.AMINO_ACIDS:
			return False
	return True


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
		elif i[0] != "'":
			seq += i
	db, identifier, entry_name, protein_name, dic = parse_fasta_header(header)
	return header, protein_name, dic['OS'], seq
	

def find_uni_start(protein_name):
	'''
	Use regular expression to find the starting url of a protein search

	Input:
		protein_name: name of the protein

	Output:
		The url of the first page of a protein search
	'''
	d = re.sub(r'[^"[A-Za-z0-9]+',' ',protein_name)
	processed = d.replace(" ","+")
	return 'https://www.uniprot.org/uniprot/?query=' + processed + \
	'&sort=score'


def find_nextpage(url):
	'''
	Takes an URL and find the url for next page

	Inputs:
		url: the current url

	Return: the url for next page
	'''
	pm = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs =\
	 certifi.where())
	html = pm.urlopen(url=url, method='GET').data
	soup = bs4.BeautifulSoup(html, "html.parser")
	tags = soup.find_all("a", class_="nextPageLink")
	if tags:
		link = set()
		for tag in tags:
			if tag.has_attr('href'):
				rv = pa2util.convert_if_relative_url(url, tag['href'])
				link.add(rv)
		return list(link)[0]
	return None

				
def get_similar(protein_name, nmax = 20):
	'''
	Find a customized number of similar proteins on the UniProt website

	Input:
		protein_name: name of the protein
		nmax: the maximum number of similar proteins that we want

	Return:
		A list of proteins UniPort codes of similar proteins 
	'''
	protein_name = protein_name.lower()
	url = find_uni_start(protein_name)
	n = 0
	result = []
	while n < nmax:
		if url == None:
			break
		similar = code_search(url)
		result += similar
		n += len(similar)
		url = find_nextpage(url)
	if n < nmax:
		diff = nmax-diff
		result += code_search(url)[:diff]		
	return result[0: nmax]


def create_MSA(similars, max_len):
	'''
	Create Multiple Sequence Alignment from the similar proteins scrapped from 
	the UniProt website and save it as fasta files.

	Input:
		similars: A list of proteins UniPort codes of similar proteins

	Return:
		None
	'''
	alignment_list = {}
	record_list = []
	for id_ in similars:
		fasta = get_fasta(id_)
		header, name, species, sequence = parse_fasta(fasta)
		db, identifier, entry_name, protein_name, dic = \
		 parse_fasta_header(header)
		if len(sequence) < max_len:
			print(len(sequence), max_len)
			record = SeqRecord(Seq(sequence, IUPAC.protein), id = identifier, 
				name = name, description = header)
			record_list.append(record)
	with open('./static/unaligned.fasta', 'w') as output_handle:
		SeqIO.write(record_list, output_handle, 'fasta')
	cline = coc(infile = './static/unaligned.fasta', 
		outfile = './static/aligned.fasta', 
		verbose = True, auto = True, force = True)
	os.system(str(cline))
	return None


class Protein:
	def __init__(self, pdb_id, uniprot_id):
		self.pdb_id = pdb_id
		self.uniprot_id = uniprot_id
		self.fasta = get_fasta(uniprot_id)
		name, species, sequence = self.parse_fasta()
		self.name = name
		self.species = species
		self.sequence = sequence
		self.length = len(self.sequence)


	def parse_fasta(self):
		'''
		Parse the FASTA file and determine the name, species, and sequence 
		of the protein.

		Output:
			name: string
			species: string
			sequence: string
		'''
		header, name, species, sequence = parse_fasta(self.fasta)
		return name, species, sequence


	def __repr__(self):
		st = self.name + " from " + self.species + " is a " + \
		 str(self.length) \
		 + " peptides long protein with UniProt ID: " + self.uniprot_id
		return st     

