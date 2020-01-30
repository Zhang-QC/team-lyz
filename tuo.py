import bs4
import requests

def get_fasta(uniprot_id):
    uniprot_url = "https://www.uniprot.org/uniprot/" + uniprot_id + ".fasta"
    r = requests.get(uniprot_url)
    text = r.text.encode()
    soup = bs4.BeautifulSoup(text, "html5lib")
    tags = soup.find_all("body")

    return tags[0].text




