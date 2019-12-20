import os
import re
import copy
import json
from xml.etree import ElementTree
import requests


def print_pretty(aresite):
    try:
        mots = aresite["exact_motifs"]
        starts = aresite["motif_starts"]
        ends = aresite["motif_ends"]
        chrs = aresite["chromosomes"]
        strands = aresite["strands"]
        transcripts = aresite["transcripts"]
        genes = aresite["genes"]
        evh = aresite["hur_evidence"]
        evt = aresite["ttp_evidence"]
        eva = aresite["auf_evidence"]
        anno = aresite["annotation"]
        meanaccs = ['{:.4e}'.format(x) for x in aresite["meanacshort"]]
        meanaccl = ['{:.4e}'.format(x) for x in aresite["meanaclong"]]

        aresite = zip(chrs, starts, ends, mots, anno, strands, genes, transcripts, meanaccs, meanaccl, evh, evt, eva)
        # exclude other than those in 3'-UTR
        aresite = [are for are in aresite if ('Exon^3UTR' in are[4])]

        def getKey(item):
            return item[1]

        aresite = sorted(aresite, key=getKey)

        for site in aresite:
            print("\t".join(site))
    except KeyError:
        pass


def process_raw_aresite_output(filename):
    with open(filename, 'rt') as f:
        try:
            data = json.load(f)
        except json.decoder.JSONDecodeError:
            return []
        # print(data)
        # print(type(data))  # <class 'dict'>
        print(filename)
        print_pretty(data)
        intervals = get_pretty(data)
        ares_with_snps = []
        for interval in intervals:
            #print(interval)
            snps = search_for_snp(interval)
            if snps:
                print(interval)
                for snp in snps:
                    print(snp)
                are_with_snp = {
                    'are': list(interval),
                    'snps': snps
                }
                ares_with_snps.append(are_with_snp)
        return ares_with_snps


def get_pretty(aresite):
    try:
        mots = aresite["exact_motifs"]
        starts = aresite["motif_starts"]
        ends = aresite["motif_ends"]
        chrs = aresite["chromosomes"]
        strands = aresite["strands"]
        transcripts = aresite["transcripts"]
        genes = aresite["genes"]
        evh = aresite["hur_evidence"]
        evt = aresite["ttp_evidence"]
        eva = aresite["auf_evidence"]
        anno = aresite["annotation"]
        meanaccs = ['{:.4e}'.format(x) for x in aresite["meanacshort"]]
        meanaccl = ['{:.4e}'.format(x) for x in aresite["meanaclong"]]

        aresite = zip(chrs, starts, ends, mots, anno, strands, genes, transcripts, meanaccs, meanaccl, evh, evt, eva)
        # exclude other than those in 3'-UTR
        aresite = [are for are in aresite if ('Exon^3UTR' in are[4])]
        return aresite
    except KeyError:
        return []


def search_for_snp(interval):
    chr = interval[0][3:]
    start = interval[1]
    end = interval[2]
    term = '"Homo sapiens"[Organism]' \
            ' AND {0}[CHR] AND {1}:{2}[CHRPOS]'.format(chr, start, end)
    URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/' \
      'esearch.fcgi'
    params = {
        'db': 'snp',
        'term': term
    }
    res = requests.get(URL, params=params)

    snps_found = []
    root = ElementTree.fromstring(res.text)
    if root.find('Count').text != '0':
        for id_element in root.iter('Id'):
            link = 'https://www.ncbi.nlm.nih.gov/snp/' + id_element.text
            snps_found.append(link)
    return snps_found


def _main():
    #term = '"Homo sapiens"[Organism] AND (C2H2[Domain Name] OR C2HC[Domain Name] OR CCHC[Domain Name])'
    term = '"Homo sapiens"[Organism] AND (C2HC[Domain Name] OR CCHC[Domain Name])'

    dirname = re.sub(r'[\\/:*?"><|]', '_', term)
    os.chdir('./' + dirname)

    with open('genes_dump.json', 'rt') as f:
        ZFPs_json = json.load(f)         # for detailed output version
    ZFPs_json_simple = copy.deepcopy(ZFPs_json)  # for simplified output version

    for i in range(len(ZFPs_json[:])):  # you can add [from:to]
        number = ZFPs_json[i]['a) number']
        ensembl = ZFPs_json[i]['c) ensembl']
        symbol = ZFPs_json[i]['b) symbol']
        filename = '{:_<4}{}_{}.json'.format(number, ensembl, symbol)

        ares_with_snps = process_raw_aresite_output(filename)
        if ares_with_snps:
            ZFPs_json[i]['e) ares_with_snps'] = ares_with_snps # add new key: value

            snps_in_ares = []
            for are in ares_with_snps:
                snps_in_ares.extend(are['snps'])
            ZFPs_json_simple[i]['e) snps_in_ares'] = snps_in_ares  # add new key: value

    with open('genes_dump_with_snps.json', 'wt') as f:
        json.dump(ZFPs_json, f, indent=4, sort_keys=True)  # detailed output
    with open('genes_dump_with_snps(simplified).json', 'wt') as f:
        json.dump(ZFPs_json_simple, f, indent=4, sort_keys=True)  # simplified output
    # result - pretty json


if __name__ == '__main__':
    _main()