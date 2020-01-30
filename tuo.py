import bs4
import requests

def get_fasta(uniprot_id):
    uniprot_url = "https://www.uniprot.org/uniprot/" + uniprot_id + ".fasta"
    r = requests.get(uniprot_url)
    text = r.text.encode()
    return text

