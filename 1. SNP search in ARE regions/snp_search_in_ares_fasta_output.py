"""This script performs the following:

1) extracts data about ARE intervals (on chr 8) from file
   (overlapping intervals are merged);
2) checks for SNPs in each interval;
if any SNPs are found in the interval:
   * their UIDs are saved - along with respective interval coords;
   * link is formed, printed (1 link per interval)
     and saved;
3) stores search results on the Entrez history server
   (therefore later use is possible).
   The subsequent search results are being associated
   with existing Search Results and saved in the same WebEnv
   (Each request is made with the same WebEnv parameter)
4) downloads fasta files for each SNP that was found,
   using data saved on Entrez history server

"""
import requests
import re
import os
from xml.etree import ElementTree
from merge import merge_overlapping  # from current directory


are_coord = []
with open('klf10_are_coord_list.txt', 'rt') as f:
    for line in f:
        are_coord.append(tuple(map(int, line.split())))  # start, end
print(are_coord)  #
are_coord = merge_overlapping(are_coord)
# merging allows:
# 1) to do less job
# 2) pick each UID only once
print(are_coord, '(merged intervals)')  #

URL_1 = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/' \
        'esearch.fcgi'
params = {
    'db': 'snp',
    'term': 'KLF10[Gene] AND "Homo sapiens"[Organism]' \
            ' AND 8[CHR] AND {start}:{end}[CHRPOS]',
    'usehistory': 'y'
}
id_found = []
links_found = []

for i in range(len(are_coord)):
    curr_params = params.copy()  # curr_params = params    won't work!
    curr_params['term'] = params['term'].format(
        start=str(are_coord[i][0]), end=str(are_coord[i][1]))
    print('# {:<4}'.format(i + 1), end='')
    res = requests.get(URL_1, params=curr_params)
    # print(res.status_code)  # 200 in all cases

    root = ElementTree.fromstring(res.text)
    print('Count: {:<5}'.format(root.find('Count').text), end='')
    print(curr_params['term'])  # Entrez query

    if root.find('Count').text != '0':
        query_key = root.find('QueryKey').text  # or str(i + 1)
        for id_element in root.iter('Id'):
            id_found.append((id_element.text, are_coord[i], query_key))  # append a tuple
        # to maintain connection between the found SNP
        # and the respective ARE and QueryKey

        # create link to html instead of xml
        link = 'https://www.ncbi.nlm.nih.gov/snp?' + res.url.partition('?')[2]
        # stage 1 - change domain
        link = re.sub(r'&?(WebEnv[^&]*|db=[^&]*|usehistory=y)&?', '', link)
        # stage 2 - remove unnecessary params
        print(link)  # link for humans is ready
        links_found.append(link)

    if i == 0:
        WebEnv = root.find('WebEnv').text
        params['WebEnv'] = WebEnv
        # adding new parameter to be used in subsequent searches

print(id_found)
print('UIDs found:')
for id_ in id_found:
    print('{:<12}  (ARE coordinates = {})'.format(id_[0], id_[1]))

URL_2 = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/' \
        'efetch.fcgi'

if not os.path.exists('./SNPs_found'):
    os.makedirs('./SNPs_found')
os.chdir('./SNPs_found')

for i in range(len(id_found)):
    params = {
        'db': 'snp',
        'WebEnv': WebEnv,
        'query_key': id_found[i][2],
        'rettype': 'fasta',
        'retmode': 'text'
    }
    print('Downloading file {} of {}'.format(
        i + 1, len(id_found)))
    res = requests.get(URL_2, params=params)
    filename = str('SNP_rs' + id_found[i][0] + '.fasta')
    with open(filename, 'wt') as f:
        f.write(res.text)

