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
    link = set()
    for tag in tags:
        if tag.has_attr('href'):
            link.add(tag['href'])
    return list(link)[0]
  




