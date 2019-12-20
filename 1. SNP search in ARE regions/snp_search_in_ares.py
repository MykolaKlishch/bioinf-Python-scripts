"""
This script performs the following:
1) extracts data about ARE intervals (on chr 8) from file;
2) checks for SNPs in each interval;
if any SNPs are found in the interval:
   * their UIDs are saved - along with respective interval coords;
   * link is formed, printed (1 link per interval)
     and saved
3) generates XML files with results (seems to be odd -
   useful only for debug. Currently behind #)
4) outputs saved data for SNPs

The program stores search results on the Entrez history server
(therefore later use is possible)
The subsequent search results are being associated
with existing Search Results and saved in the same WebEnv
(Each request is made with the same WebEnv parameter)

But:
* since SNP UID can be used to download only XML about this SNP
  (not sth more interesting like .gb or .fasta file)
  obtained UIDs are just printed and not used anymore
* WebEnv is not used anymore
* in this case links are more useful output than UIDs or XMLs
History usage was added just for training while writing this script.
This option is useless in the present script
But it can be useful if this script is used to create another script
in which UIDs history really matters.

Currently history is disabled to shorten links
"""

import requests
from xml.etree import ElementTree
# import time


are_coord = []
with open('klf10_are_coord_list.txt', 'rt') as f:
    for line in f:
        are_coord.append(':'.join(line.split()))  # start, end
print(are_coord)

URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/' \
      'esearch.fcgi'
params = {
    'db': 'snp',
    'term': 'KLF10[Gene] AND "Homo sapiens"[Organism]' \
            ' AND 8[CHR] AND {}[CHRPOS]',
    # 'usehistory': 'y'
}
id_found = []
links_found = []

for i in range(len(are_coord)):
    curr_params = params.copy()  # curr_params = params    won't work!
    curr_params['term'] = params['term'].format(are_coord[i])
    print('# {:<4}'.format(i + 1), end='')
    res = requests.get(URL, params=curr_params)
    # print(res.status_code)  # 200 in all cases
    print('url: ', res.url)  # not interested in this url

    # filename = 'snp_in_are{}_usehistory=y'.format(str(i + 1))
    # with open(filename + '.xml', 'w') as f:
    #     f.write(res.text)                     # -- USELESS !!!

    root = ElementTree.fromstring(res.text)
    print('Count: {:<5}'.format(root.find('Count').text), end='')
    print(curr_params['term'])  # Entrez query

    if root.find('Count').text != '0':
        for id_element in root.iter('Id'):
            id_found.append((id_element.text, are_coord[i]))  # append a tuple
        # to maintain connection between the found SNP and the respective ARE
        link = 'https://www.ncbi.nlm.nih.gov/snp?' + res.url.partition('?')[2]
        print(link)
        links_found.append(link)
        # print('***** SNP found! *****')
    # if i == 0:
    #     WebEnv = root.find('WebEnv').text
    #     params['WebEnv'] = WebEnv
    #   # adding new parameter to be used in subsequent searches
    # time.sleep(0.35)  # in order not to be banned - just in case

print('UIDs found:')
for id_ in id_found:
    print('{:<12}  (ARE coordinates = {})'.format(id_[0], id_[1]))
# print('Links o all found SNPs')
# for link in links_found:  -- USELESS !!!
#     print(link)
