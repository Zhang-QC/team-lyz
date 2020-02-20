import bs4
import requests
import re
import urllib3
import certifi
import pa2util



AMINO_ACIDS = {"A": "alanine", "R": "arginine", "N": "asparagine", 
"D": "aspartate", "C": "cysteine", "E": "glutamate", "Q": "glutamine",
"G": "glycine", "H": "histidine", "I": "isoleucine", "L": "leucine",
"K": "lysine", "M": "methionine", "F": "phenylalanine", "P": "proline",
"S": "serine", "T": "threonine", "W": "tryptophan", "Y": "tyrosine",
"V": "valine"}

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

def get_fasta(uniprot_id):
    uniprot_url = "https://www.uniprot.org/uniprot/" + uniprot_id + ".fasta"
    r = requests.get(uniprot_url)
    text = r.text.encode()
    soup = bs4.BeautifulSoup(text, "html5lib")
    tags = soup.find_all("body")

    return str(tags[0].text)


def read_fasta(fast_str):
    strin = fast_str.split('\n')
    substr = strin[0]
    sub_list = re.find_all('([A-Z][0-9]{5}),([A-Z0-9]{4}_[A-Z]{5}),(OS=[A-Z][a-z\s]+),(OX=[0-9]{4}),(GN=[A-Za-z0-9]{5}),(PE=[0-9]),(SV=[0-9])', substr)
    total_list = sub_list
    for sub_list in strin[1:]:
        if sub_list != ' ':
            total_list.append(sub_list)
    return total_list
    

def code_search(url):
    r = requests.get(url)
    text = r.text.encode()
    soup = bs4.BeautifulSoup(text,"html5lib")
    tags = soup.find_all("a")
    s = ""
    for tag in tags:
        s += tag.text
    match=re.findall('[A-Z][0-9][A-Z0-9]{3}[0-9]',s)
    return match



def find_uni_start(protein_name):
    return 'https://www.uniprot.org/uniprot/?query=' + protein_name +'&sort=score'

def find_nextpage(url):
    '''
    Takes an URL and find the url for next page

    Inputs:
        url: the current url

    Return: the url for next page
    '''
    pm = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs = certifi.where())
    html = pm.urlopen(url=url, method='GET').data
    soup = bs4.BeautifulSoup(html)
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
    protein_name = protein_name.lower()
    url = find_uni_start(protein_name)
    n = 0
    result = []
    while n < nmax:
        if url == None:
            break
        similar = code_search(url)
        result += similar
        n+=len(similar)
        url = find_nextpage(url)
        print(url)
    if n < nmax:
    	diff = nmax-diff
    	result+= code_search(url)[:diff]
    return result



class protein:
	def __init__(self, pdb_id, uniprot_id):
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

		Output:
			name: string
			species: string
			sequence: string
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

	def __repr__(self):
		str = self.name + " from " + self.species + " is a " + self.length \
		+ "peptides long protein with UniProt ID: " + self.uniprot_id
		return str		









































































