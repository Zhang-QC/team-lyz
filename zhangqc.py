'''
***********
Jason Zhang
***********
'''

import bs4
import requests
import util
import tuo
import jenny


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


