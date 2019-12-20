import os
import json
import re
import time
from final_part_2 import get_aresite_json


def _main():
    term = '"Homo sapiens"[Organism] AND (C2H2[Domain Name] OR C2HC[Domain Name] OR CCHC[Domain Name])'
    # term = '"Homo sapiens"[Organism] AND (C2HC[Domain Name] OR CCHC[Domain Name])'

    dirname = re.sub(r'[\\/:*?"><|]', '_', term)
    os.chdir('./' + dirname)

    with open('genes_dump.json', 'rt') as f:
        znfs_json = json.load(f)

    invalid_files = []
    for znf in znfs_json:
        number = znf['a) number']
        ensembl = znf['c) ensembl']
        symbol = znf['b) symbol']
        filename = '{:_<4}{}_{}.json'.format(number, ensembl, symbol)
        print(filename)

        reason = ''
        with open(filename, 'rt') as f:
            file_content = f.read()
            #print(file_content)
            if file_content == '':
                reason = 'empty file'
            if '<!DOCTYPE html>' in file_content:
                reason = 'file contains error report'
            # if '"reason":"access denied"' in file_content:
            #     reason = 'access denied'
        if reason:
            print('    #####    ', reason)
            znf['reason'] = reason
            invalid_files.append(znf)

    with open('invalid_files.json', 'wt') as f:
        json.dump(invalid_files, f, indent=4, sort_keys=True)

    print('\n{} INVALID FILES FOUND\n'.format(str(len(invalid_files))))

    timeouts_exceeded = 0  # in a row

    # invalid_files_part = invalid_files[:]
    start = int(input('start from number: '))
    for znf in invalid_files:
        print(znf['a) number'], end=' ')
        if int(znf['a) number']) >= start:
            print('\n', invalid_files.index(znf))
            invalid_files_part = invalid_files[invalid_files.index(znf):]
            break

    for znf in invalid_files_part:
        number = znf['a) number']
        ensembl = znf['c) ensembl']
        symbol = znf['b) symbol']
        filename = '{:_<4}{}_{}.json'.format(number, ensembl, symbol)
        print(filename)
        try:
            start_time = time.clock()

            raw_json = get_aresite_json(ensembl, 'Homo_sapiens', 'ATTTA')
            with open(filename, 'wt') as f:
                f.write(raw_json)

            timedelta = time.clock() - start_time
            if timedelta >= 120:
                timeouts_exceeded += 1
            else:
                timeouts_exceeded = 0
            print("request processed in {:<2.4} seconds".format(str(timedelta)))

            time.sleep(timedelta)  # To avoid overload
            # The slower is the server, the more time it has to rest
        except Exception as e:
            print(e.args)
            print(e)
            print('Could not properly dump ARE info')
        if timeouts_exceeded >= 3:  # 3 in a row
            break

if __name__ == "__main__":
    _main()