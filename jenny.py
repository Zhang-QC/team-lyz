'''
**********
Jenny Yang
**********
'''
import bs4

import requests
import re
import util
import bs4
import queue
import json
import sys
import csv
import pa2util
import urllib3
import certifi

from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment



def find_nextpage(url):
    '''
    Takes an URL and find the url for next page

    Inputs:
        url: the current url

    Return: the url for next page
    '''
    pm = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs = certifi.where())
    html = pm.urlopen(url=url, method='GET').data
    soup = bs4.BeautifulSoup(html,features="lxml")
    tags = soup.find_all("a", class_="nextPageLink")
    link = set()
    for tag in tags:
        if tag.has_attr('href'):
            link.add(tag['href'])
    return list(link)[0]


def find_uni_start(protein_name):
    return 'https://www.uniprot.org/uniprot/?query=' + protein_name +'&sort=score'

#prefer to have dictionary with id as key and alignment as value
def create_MSA(alignment_list):
    align = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
    for id_ in alignment_list.keys():
         align.add_sequence(id_, alignment_list[id_])
    print(align)



  




