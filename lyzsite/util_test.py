import bs4
import requests
import re
import urllib3
import certifi
import pa2util
import main



des = {"first": "Percent identity of seqences within the multiple sequence \
alignment",
"2": "Conservation of residues in the reference sequence",
"3": "Correlation matrix of the Statistical Coupling Analysis",
"4": "Conservation map of protein onto the 3D protein model (Front)",
"5": "Conservation map of protein onto the 3D protein model (Back)",
}


def get_fasta(args):
    '''
    Take django input and run our program

    Inputs:
        agrs: dictionary created by Django using users' input

    Returns: Pictures and descriotions
    '''
    if not args:
        return [[],[]]
    name = args["terms"]
    length = args['Maximum number of sequences']
    evalue = args['L value']
    main.run_all(name, length, evalue)
    return [['picture','des'], [("1.png", des['first']),("2.png", des['2']),\
    ("3.png", des['3']),("4.png", des['4']), ("5.png", des['5'])]]












































