import bs4
import requests

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


class protein:
	def __init__(self, pdb_id = '', uniprot_id):
		self.pdb_id = pdb_id
		self.uniprot_id = uniprot_id
		self.fasta = tuo.get_fasta(uniprot_id)
		name, species, sequence = self.parse_fasta()
		self.name = name
		self.species = species
		self.sequence = sequence
		self.length = len(self.sequence)

	def parse_fasta(self):
		'''
		Parse the FASTA file and determine the name, species, and sequence 
		of the protein.
		'''
		name = ''
		species = ''
		sequence = ''
		return name, species, sequence

	def find_similar(self, len_diff, max_num, curated = True):
		'''
		Find similar proteins from the UniProt database using the assigned 
		parameters.
		'''
		return None

	def blast(self):
		'''
		Perform a BLAST search on the protein sequence.
		'''
		seq = self.sequence
		return None