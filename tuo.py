import bs4
import requests

def get_fasta(uniprot_id):
    uniprot_url = "https://www.uniprot.org/uniprot/" + uniprot_id + ".fasta"
    r = requests.get(uniprot_url)
    text = r.text.encode()
    soup = bs4.BeautifulSoup(text, "html5lib")
    tags = soup.find_all("body")

    return tags[0].text




def read_fasta(fast_str):
    strin = fast_str.split('\n')
    substr = strin[0]
    sub_list = re.find_all('([A-Z][0-9]{5}),([A-Z0-9]{4}_[A-Z]{5}),(OS=[A-Z][a-z\s]+),(OX=[0-9]{4}),(GN=[A-Za-z0-9]{5}),(PE=[0-9]),(SV=[0-9])', substr)
    total_list = sub_list
    for sub_list in strin[1:]:
        if sub_list != ' '
            total_list.append(sub_list)
    return total_list
