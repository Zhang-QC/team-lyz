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

def find_next(url):
    next_page_text = bs.find('ulr', class_="SearchBreadcrumbs").findAll('li')[-1].text
    if next_page_text = 'Next':
        next_page_partial = bs.find('ul',class_="SearchBreadcrumbs").findAll('li')[-1].find('a')['href']
            next_page_url = base_url + next_page_partial
            return(next_page_url)
    else:
        return None


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

