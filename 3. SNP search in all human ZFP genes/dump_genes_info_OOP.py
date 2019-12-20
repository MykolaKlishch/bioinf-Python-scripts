import os
import re
import sys
import json
from xml.etree import ElementTree
import requests


class Gene:

    def __init__(self, raw_info):
        self.raw_info = raw_info
            # raw info about a single gene
            # do we need this attr - ?
        match_obj = re.match(r'\b([\d]+)\. (.+)\b', raw_info)
        self.number = match_obj.group(1)
        self.symbol = match_obj.group(2)
        self.ensembl = get_ensembl(self.symbol)
            # if the gene with such symbol cannot be found in ensembl db
            # self.ensembl will be None
            # pay attention on this when using self.ensembl attribute!
        self.chr = re.search(r'Chromosome: (\w+)', raw_info).group(1)
            # not (\d+) - chr is not always a number - it can also be X or Y
        # try:
        #     self.number = match_obj.group(1)
        # except AttributeError:
        #     self.number = 'Not found'
        # try:
        #     self.symbol = match_obj.group(2)
        # except AttributeError:
        #     self.symbol = 'Not found'
        # if self.symbol:
        #     self.ensembl = self.get_ensembl(self.symbol)
        # else:
        #     self.ensembl = 'Impossible to obtain'
        # try:
        #     are_intervals = get_intervals(
        #     self.ensembl, 'Homo_sapiens', 'ATTTA'
        #     # Homo_sapiens works. But homo_sapiens - doesn't!!!
        # )
        # except KeyError:
        #     are_intervals = ['Error occured']
        # self.are_intervals = are_intervals


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
        #'term': '{}[uid]'.format(self.uid),
        'WebEnv': WebEnv,
        'query_key': '1',
        'rettype': 'json',
        'retmode': 'text'
    }
    res = requests.get(URL, params=params)
    return res.text.strip()  # strip is important for correct further split

def dump_genes_info(info, term):
    """Creates new directory; its name corresponds to query.
    Outputs a list of gene UIDs as a txt file in this directory

    :param info: text
    :return: None

    """
    dirname = re.sub(r'[\\/:*?"><|]', '_', term)
    filename = 'dump_info_in_txt.txt'

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


def get_intervals(gene, species='Homo_sapiens', motifs='ATTTA'):
    URL = 'http://rna.tbi.univie.ac.at/AREsite2/api/'
    params = {
        'query': gene,
        'species': species,
        'list': motifs
    }
    headers = {
        'Content-Type': 'application/json'
    }
    try:
        res = requests.get(URL, params=params, headers=headers, timeout=60)
    except requests.exceptions.ReadTimeout:
        return []
    print(res.text)
    print('))))))))))))))))))))))))))))))))))))))))')
    data = res.json()
    print(data)
    if data:
        #print(data.keys())
        if "motif_starts" in data.keys() and "motif_ends" in data.keys():
            starts = [int(start) for start in data["motif_starts"]]
            ends = [int(end) for end in data["motif_ends"]]
            intervals = list(zip(starts, ends))
        else:
            intervals = []
        return intervals


def get_aresite_json(gene, species='Homo_sapiens', motifs='ATTTA'):
    URL = 'http://rna.tbi.univie.ac.at/AREsite2/api/'
    params = {
        'query': gene,
        'species': species,
        'list': motifs
    }
    headers = {
        'Content-Type': 'application/json'
    }
    try:
        res = requests.get(URL, params=params, headers=headers, timeout=60)
    except requests.exceptions.ReadTimeout:
        return '[]'
    return (res.text)


def _main():
    # are_intervals = get_intervals('ZNF474')
    # for interval in are_intervals:
    #     print(interval)
    term = '"Homo sapiens"[Organism] AND (C2H2[Domain Name] OR C2HC[Domain Name] OR CCHC[Domain Name])'
    #term = '"Homo sapiens"[Organism] AND (C2HC[Domain Name] OR CCHC[Domain Name])'
    WebEnv = get_WebEnv(term)
    raw_info = get_genes_info(WebEnv) # raw info about all the genes - as single text
    dump_genes_info(raw_info, term)
    raw_info_list = [notation.strip() for notation in raw_info.split('\n\n')] # raw info about all the genes - as list
    # stripping the elements is necessary
    # because some of them are separated not by \n\n but by \n\n\n
    # in this case splitting may result in elements starting with (extra) \n
    # extracting number and symbol from such piece of raw data would raise error

    ZNFs = []  # list of Gene objects
    for descr in raw_info_list[42:]:
        print(descr)
        ZNF = Gene(descr)
        ZNFs.append(ZNF)
        # print(next_ZNF.number)
        # print(next_ZNF.chr)
        # print(next_ZNF.symbol)
        print(ZNF.ensembl)
        try:
            raw_json = get_aresite_json(ZNF.ensembl, 'Homo_sapiens', 'ATTTA')
            filename = '{:_<4}{}_{}.json'.format(ZNF.number, ZNF.ensembl, ZNF.symbol)
            with open(filename, 'wt') as f:
                f.write(raw_json)
            # with open(filename, 'wt') as f:
            #     json.dump(raw_json.json(), f, indent=4, sort_keys=True)
            # # with open(filename, 'rt') as f:
            # #     print(json.load(f))
        except Exception:
            print('Could not properly dump ARE info')
        print('++++++++')

    print(ZNFs)
    print(raw_info_list)


if __name__ == '__main__':
    _main()