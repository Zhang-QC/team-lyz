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

def find_next(url):
    next_page_text = bs4.find('url', class_="SearchBreadcrumbs").findAll('li')[-1].text
    if next_page_text == 'Next':
        next_page_partial = bs4.next_sibling('ul',class_="SearchBreadcrumbs").findAll('li')[-1].next_sibling('a')['href']
        next_page_url = base_url + next_page_partial
        return(next_page_url)
    else:
        return None
def find_nextpage(url):
    request = pa2util.get_request(url)
    links = []
    soup = None
    
    if request:
        re_url = pa2util.get_request_url(request)
    
    
    soup = bs4.BeautifulSoup(pa2util.read_request(request),'html5lib')
    tags = soup.find_all("a", class_="nextPageLink")
    link = set()
    for tag in tags:
        if tag.has_attr('href'):
            link.add(tag['href'])
    return link
  







def go_next(max_page,result_filename):
    starting_url = ('')
    num_visit = 0
    match = []
    url = starting_url
    while num_visit < max_page and final:
        url = find_next(url)
        if url:
            match.append(code_search(url))
            num_visit += 1
        else:
            final = True
    return match

