import os
import re
import sys
import json
from xml.etree import ElementTree
import requests


class Gene:

    def __init__(self, raw_info):
        self.raw_info = raw_info
        match_obj = re.match(r'\b([\d]+)\. (.+)\b', raw_info)
        self.number = match_obj.group(1)
        self.symbol = match_obj.group(2)
        self.ensembl = get_ensembl(self.symbol)
            # if the gene with such symbol cannot be found in ensembl db
            # self.ensembl will be None
            # pay attention on this when using self.ensembl attribute!
        self.chr = re.search(r'Chromosome: (\w+)', raw_info).group(1)
            # not (\d+) - chr is not always a number - it can also be X or Y


def get_WebEnv(term, retmax=1000):
    """Performs the search using history
    Stores gene UIDs on the server
    Returns WebEnv value for subsequent usage

    :param term: string
        Query (in NCBI advanced search syntax)
    :param retmax: int
        Maximum number of items in the result.
        Others are being truncated.
    :return: WebEnv for subsequent usage in efetch

    """
    URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/' \
          'esearch.fcgi'
    params = {
        'db': 'gene',
        'term': term,
        'retmax': str(retmax), #default value for NCBI search is only 20!
        'usehistory': 'y',
        'sort': 'name'
    }
    res = requests.get(URL, params=params)
    root = ElementTree.fromstring(res.text)
    WebEnv = root.find('WebEnv').text
    return WebEnv


def get_genes_info(WebEnv):
    """
    :param WebEnv: web environment for fetching
    :return: Returns text containing information about all the genes
    corresponding to the term

    """
    URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/' \
          'efetch.fcgi'
    params = {
        'db': 'gene',
        'WebEnv': WebEnv,
        'query_key': '1',
        'rettype': 'json',
        'retmode': 'text'
    }
    res = requests.get(URL, params=params)
    return res.text.strip()  # strip is important for correct further split


def dump_genes_info(info, term):
    """Creates new directory; its name corresponds to query.
    Outputs a text containing information about all the genes
    as a txt file in this directory

    :param info: text
    :return: None

    """
    dirname = re.sub(r'[\\/:*?"><|]', '_', term)
    filename = 'genes_dump.txt'

    if not os.path.exists('./' + dirname):
        os.makedirs('./' + dirname)
    os.chdir('./' + dirname)

    with open(filename, 'wt') as f:
        f.write(info)


def get_ensembl(symbol):
    """Sends request containing gene symbol (e.g. KLF4)
    to the server. Receives ENSEMBL code

    """
    URL = "http://rest.ensembl.org/" \
          "xrefs/symbol/homo_sapiens/{}?".format(symbol)
    headers = {
        "Content-Type": "application/json"
    } # params won't work
    res = requests.get(URL, headers=headers)

    if not res.ok:
        res.raise_for_status()
        sys.exit()

    decoded = res.json()
    for obj in decoded:
        if obj['id'].startswith('ENS'):
            return obj['id']


def _main():
    #term = '"Homo sapiens"[Organism] AND (C2H2[Domain Name] OR C2HC[Domain Name] OR CCHC[Domain Name])'
    term = '"Homo sapiens"[Organism] AND (C2HC[Domain Name] OR CCHC[Domain Name])'
    WebEnv = get_WebEnv(term)
    raw_info = get_genes_info(WebEnv)  # raw info about all the genes - as single text
    dump_genes_info(raw_info, term)  # 1. dump raw info in txt
    raw_info_list = [descr.strip() for descr in raw_info.split('\n\n')] # raw info about all the genes - as list
    # stripping the elements is necessary
    # because some of them are separated not by \n\n but by \n\n\n
    # in this case splitting may result in elements starting with (extra) \n
    # extracting number and symbol from such piece of raw data would raise error

    ZNFs = []  # list of Gene objects
    ZNFs_json = []  # list of dicts to become json
    for descr in raw_info_list:
        ZNF = Gene(descr)
        ZNFs.append(ZNF)
        print(ZNF.raw_info)
        print(ZNF.ensembl)
        ZNF_dict = {
            'a) number': ZNF.number,
            'b) symbol': ZNF.symbol,
            'c) ensembl': ZNF.ensembl,
            'd) chr': ZNF.chr
        }
        ZNFs_json.append(ZNF_dict)

    with open('genes_dump.json', 'wt') as f:
        json.dump(ZNFs_json, f, indent=4, sort_keys=True)
        # 2. dump info in json


if __name__ == '__main__':
    _main()