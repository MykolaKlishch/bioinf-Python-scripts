import os
import re
import json
import time
import requests


def get_aresite_json(gene, species='Homo_sapiens', motifs='ATTTA'):
    url = 'http://rna.tbi.univie.ac.at/AREsite2/api/'
    params = {
        'query': gene,
        'species': species,
        'list': motifs
    }
    headers = {
        'Content-Type': 'application/json'
    }
    try:
        res = requests.get(url, params=params, headers=headers, timeout=360)  # ####
    except requests.exceptions.ReadTimeout:
        return ''
    return res.text


def _main():
    term = '"Homo sapiens"[Organism] AND (C2H2[Domain Name] OR C2HC[Domain Name] OR CCHC[Domain Name])'
    # term = '"Homo sapiens"[Organism] AND (C2HC[Domain Name] OR CCHC[Domain Name])'

    dirname = re.sub(r'[\\/:*?"><|]', '_', term)
    os.chdir('./' + dirname)

    with open('genes_dump.json', 'rt') as f:
        znfs_json = json.load(f)

    start = int(input('Start from number: '))
    timeouts_exceeded = 0  # in a row

    for i in range(start - 1, len(znfs_json), 1):
        number = znfs_json[i]['a) number']
        ensembl = znfs_json[i]['c) ensembl']
        symbol = znfs_json[i]['b) symbol']
        filename = '{:_<4}{}_{}.json'.format(number, ensembl, symbol)
        print(filename)
        try:
            start_time = time.clock()

            raw_json = get_aresite_json(ensembl, 'Homo_sapiens', 'ATTTA')
            with open(filename, 'wt') as f:
                f.write(raw_json)

            timedelta = time.clock() - start_time
            if timedelta >= 360:   # ####
                timeouts_exceeded += 1
            else:
                timeouts_exceeded = 0
            print("request processed in {:<2.4} seconds".format(str(timedelta)))

            time.sleep(timedelta)   # To avoid overload
            # The slower is the server, the more time it has to rest
        except Exception as e:
            print(e.args)
            print(e)
            print('Could not properly dump ARE info')
        if timeouts_exceeded >= 3:  # 3 in a row
            break


if __name__ == '__main__':
    _main()
