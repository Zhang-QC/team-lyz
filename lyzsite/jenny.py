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


# The following code is related to calculating SCA

q = [0.073,0.025,0.050,0.061,0.042,0.072,0.023,\
0.053,0.064,0.089, 0.023, 0.043, 0.052, 0.040,\
 0.052, 0.073, 0.056, 0.063, 0.013, 0.033]

def get_frequency(align, amino, position):
    '''
    Calculate the conservation of an amino acid at a certain position

    Input:
        align: multiple seqence alignment object
        amino: the target amino acid
        position: the freqeuncy at at certain position within an alignment

    Return: 
        conservation of an amino acid at a position
    '''
    M = len(align)
    n = 0
    for record in align:
        if record.seq[position] == amino:
            n+=1
    return n/M



def get_conservation(align, amino, position):
    '''
    Calculate the conservation of an amino acid at a certain position

    Input:
        align: multiple seqence alignment object
        amino: the target amino acid
        position: the freqeuncy at at certain position within an alignment

    Return: 
        conservation of an amino acid at a position
    '''
    M = len(align)
    f = get_frequency(align, amino, position)
    #what if q len different from MSA
    return f*ln(f/q[position])+(1-f)*ln((1-f)/(1-q[position]))

def joint_dis(align, a1, a2, p1, p2):
    M = len(align)
    n = 0
    for record in align:
        if record.seq[p1] == a1 and record.seq[p1] == a2:
            n+=1
    return n/M




def cov(align, a1, a2, p1, p2):
    f1 = get_frequency(a1,p1)
    f2 = get_frequency(a2,p2)
    f12 = joint_dis(align, a1, a2, p1, p2)
    w1 = abs(ln(f1*(1-q[p1])/(q[p1]*(1-f1))))
    w2 = abs(ln(f2*(1-q[p2])/(q[p2]*(1-f2))))
    return w1*w2*(f12-f1*f2)





