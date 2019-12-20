# Copyright (c) 2019 Mykola Klishch. All rights reserved. No warranty.

"""
This script performs the following:

1. Extracts information from GenBank cDNA file and creates FASTA format file.
   The latter one contains mRNA sequence with coding region and AREs highlited
   in uppercase. The rest of the sequence remains in lowercase.

2. Extracts information from FASTA file about SNPs found in  this region.
   To obtain this FASTA file, you should do the following:
    > go to https://www.ncbi.nlm.nih.gov/snp
    > perform SNP search in this region
      It's recommended that search region corresponds with the sequence
      in from GenBank mRNA file (see above). If this file contains SNPs
      from different region, the program won't be able to mark them in
      the mRNA (cDNA) sequence obtained from GenBank file.
    > save search results as follows:
      Send To > File > Format: FASTA > Create file

3. Marks SNPs in the mRNA (cDNA) sequence. It is done using letters for
   ambiguous nucleotides. As a result, consensus sequence is obtained.
   This sequence contains "ambiguous letters" (e.g. R, Y, W, etc.)
   in all polymorphic positions.

4. Highlights AREs (and coding region) in consensus mRNA (cDNA) sequence
   (in the same way, as in (1)).

5. Compares highlights in original and SNP-consensus sequences.
   Comparing highlights allows to find "newly emerged AREs"
   i.e. AREs that are present only in SNP minor allele carriers.
   UIDs (rs) of these SNPs are reported in the output.
"""

import re
import requests

class SNP:

    def __init__(self, snp_fasta_item):
        self.uid = None         # UID of the SNP - the same as rs
        self.flank_5 = None     # 5'-flanking sequence
        self.flank_3 = None     # 3'-flanking sequence
        self.clas = None        # class: snp or in-del
                                # indels can be found among SNPs
        self.alleles = None     # e.g. ('A', 'G')
        self.flanks_were_reversed = False
                                # in case given flanking sequences
                                # were not from coding strand
                                # it's useful to reverse them
        self.position = None    # position in this particular sequence
        self.chrpos = None      # position on the whole chromosome

        try:
            self.uid = re.search(
                r"""dbSNP\|(rs\d*?) """,
                snp_fasta_item).groups()[0].strip()
            flank_5_raw, flank_3_raw = re.search(
                r"""\n([ATGC \n]*?)\n\w{1}\n([ATGC \n]*?)\n\n""",
                snp_fasta_item).groups()
            # these sequences are considered raw
            # because they contain ' ' and '\n'
            self.flank_5 = ''.join(flank_5_raw.strip().split()).strip()
            self.flank_3 = ''.join(flank_3_raw.strip().split()).strip()
            self.flank_5 = self.flank_5 # [len(self.flank_5)-20:]######## make them shorter - as a test
            self.flank_3 = self.flank_3 # [:20]########
            # now the sequences are clear
            self.clas = re.search(
                r'''class=(.*?)\|''',
                snp_fasta_item).groups()[0].strip()
            alleles_raw = re.search(
                r'''alleles="(.*?)"''',
                snp_fasta_item).groups()[0].strip()
            self.alleles = frozenset(alleles_raw.split('/'))
        except AttributeError or ValueError:
            pass


def read_from_GenBank_file(name):
    """Reads a .gb file.

    The function extracts the following information:
    * definition of the sequence;
    * accession number (incl.version);
    * organism name;
    * protein sequence;
    * rna sequence.

    :returns
        tuple: 5 parameters (in the same order)
        which can be assigned to global values.

    """
    file = open(name + '.gb', 'rt')
    # 2 following wariables indicate whether
    # sequence reading is in progress
    prot_rd = False
    rna_rd = False
    try:
        for line in file:
            if line.startswith('DEFINITION'):
                DEF = line.partition('DEFINITION')[2].strip()
            elif line.startswith('VERSION'):
                ACC = line.partition('VERSION')[2].strip()
            elif line.startswith('  ORGANISM'):
                ORG = line.partition('ORGANISM')[2].strip()
            elif prot_rd:
                if '"' in line:  # final line of protein sequence
                    prot_seq += line.strip()[:-1]
                    prot_rd = False
                else:
                    prot_seq += line.strip()
            elif line.partition('/translation="')[2] != '':
                # /translation="" tag does not start the line
                # there are many ' ' before it
                prot_seq = line.partition('/translation="')[2].strip()
                prot_rd = True
            elif rna_rd:
                if '//' in line:  # no more lines with rna sequence
                    rna_rd = False
                else:
                    rna_seq += ''.join(line.strip().split()[1:])
            elif line.startswith('ORIGIN'):
                rna_seq = ''
                rna_rd = True
                # start reading rna sequence from the next line
        return (DEF, ACC, ORG, prot_seq, rna_seq)
    except ValueError or Exception:
        print('Unexpected file content')
    finally:
        file.close()


GENETIC_CODE = {
    'A': 'GCN',  # ala: GCT, GCC, GCA, GCG
    'R': 'MGN',  # arg: CGT, CGC, CGA, CGG, AGA, AGG
    'N': 'AAY',  # asn: AAT, AAC
    'D': 'GAY',  # asp: GAT, GAC
    'C': 'TGY',  # cys: TGT, TGC
    'Q': 'CAR',  # gln: CAA, CAG
    'E': 'GAR',  # glu: GAA, GAG
    'G': 'GGN',  # gly: GGT, GGC, GGA, GGG
    'H': 'CAY',  # his: CAT, CAC
    'I': 'ATH',  # ile: ATT, ATC, ATA
    'L': 'YTN',  # leu: TTA, TTG, CTT, CTC, CTA, CTG
    'K': 'AAR',  # lys: AAA, AAG
    'M': 'ATG',  # met: ATG (start)
    'F': 'TTY',  # phe: TTT, TTC
    'P': 'CCN',  # pro: CCT, CCC, CCA, CCG
    'S': 'WSN',  # ser: TCT, TCC, TCA, TCG, AGT, AGC
    'T': 'ACN',  # thr: ACT, ACC, ACA, ACG
    'W': 'TGG',  # trp: TGG
    'Y': 'TAY',  # tyr: TAT, TAC
    'V': 'GTN',  # val: GTT, GTC, GTA, GTG
    'X': 'NNN',  # any amino acid
}


def reverse_translate(prot_seq):
    """Returns consensus coding sequence for the protein.

    Sequences are represented in the standard IUB/IUPAC
    amino acid and nucleic acid codes, including letters
    for ambiguous nucleotides (R for A and G, N for any nucleotide etc.).
    https://en.wikipedia.org/wiki/FASTA_format#Sequence_representation

    GENETIC_CODE - <class 'dict'> constant.
    Contains inverse table for the standard genetic code (DNA version)
    Compressed using IUPAC notation.
    https://en.wikipedia.org/wiki/Genetic_code#Standard_codon_tables

    """
    rev_tr = ''
    for aminoacid in prot_seq:
        rev_tr += GENETIC_CODE[aminoacid]
    rev_tr += 'TRR'
    # consensus stop codon
    return (rev_tr)


NUCL_CODE = {
    'A': frozenset('A'),  # Adenine
    'C': frozenset('C'),  # Cytosine
    'G': frozenset('G'),  # Guanine
    'T': frozenset('T'),  # Thymine. U was removed !!!
    'R': frozenset('AG'),  # puRine
    'Y': frozenset('CT'),  # pYrimidines
    'K': frozenset('GT'),  # bases which are Ketones
    'M': frozenset('AC'),  # bases with aMino groups
    'S': frozenset('CG'),  # Strong interaction
    'W': frozenset('AT'),  # Weak interaction
    'B': frozenset('CGT'),  # not A - B comes after A
    'D': frozenset('AGT'),  # not C - D comes after C
    'H': frozenset('ACT'),  # not G - H comes after G
    'V': frozenset('ACG'),  # neither T [nor U] - V comes after U
    'N': frozenset('ACGT'),  # any Nucleotide
}


def match(nucl1, nucl2):
    """The function returns True if base inentity is NOT EXCLUDED.

    Can also compare ambiguous nucleotides in both sequences - e.g.:
        R & A or G => True
        N & any nucleotide => True
        N & R => True etc.)

    NUCL_CODE - <class 'dict'> constant.
    Contains standard IUB/IUPAC nucleic acid codes, including letters
    for ambiguous nucleotides
    https://en.wikipedia.org/wiki/FASTA_format#Sequence_representation

    NUCL_CODE resides outside the function in order to avoid
    multiple reassigning of the same value.
    match() is frequently used function in this program).

    """
    return (not NUCL_CODE[nucl1].isdisjoint(NUCL_CODE[nucl2]))


def scan_for_overlaps(seq, q_subseq, start=0):
    """Finds the region in seq that completely overlaps with cons_seq

    (i.e. match() returns True for all nucleotides).
    Returns (as tuple):
     * sequence with overlapping region in uppercase
     * index of the first nucl after overlap
     (the start of 3'-UTR if q_subseq is coding sequence).

    Function applies uppercase() to this region. Other parts are output
    in the unchanged case. Therefore input should be in lowercase
    (at least the surroundings of overlapping regions)
    in order to distinguish the overlap.

    Sequences are expected to be represented in the
    standard IUB/IUPACnucleic acid codes.
    Both sequences may include letters for ambiguous nucleotides
    https://en.wikipedia.org/wiki/FASTA_format#Sequence_representation.

    Gaps are not supported!

    """
    or_seq = seq
    # letter case in original sequense remains unchanged
    seq = seq.upper()
    q_subseq = q_subseq.upper()
    # case unification allows detection of overlapping with sequences
    # that already contain different cases. This can be useful e.g. when
    # searching for ARE motif WWTTTWW after WTTTW has already been detected
    continue_from = 0
    for i in range(start, len(seq) - len(q_subseq) + 1):
        for j in range(len(q_subseq)):
            if not match(seq[i + j], q_subseq[j]):
                break
        else:
            or_seq = or_seq[:i] + or_seq[i:i + j + 1].upper() + or_seq[i + j + 1:]
            # original sequense is changed:
            continue_from = i + j + 1
            # to determine 3'-UTR start after finding CDS.
            # Then scan for AREs from 3'-UTR start
    return (or_seq, (continue_from if or_seq != seq else 0))
    # with changed or intact letter case
    # return is used outside the cycle thus allowing multiple detection
    # of the same sequence


AREs = ['ATTTA',
        'WTTTW',
        'WWTTTWW',
        'WWWTTTWWW',
        'WWWWTTTWWWW',
        'WWWWWTTTWWWWW',
        'TTTGTTT',
        'GTTTG',
        'AWTAAA'
        ]  # source: http://rna.tbi.univie.ac.at/AREsite2/welcome


def output_as_FASTA(filename, seq, ACC, DESCR):
    """Creates .fasta file with sequence from seq
    or complements existing .fasta file.

    The file contains:
     * accession number (obligatory);
     * description(e.g. full definition or organism name);
     * nucleotide sequence.

    """
    f = open(filename + '.fasta', 'at')
    f.write('>{0} {1}\n'.format(ACC, DESCR))
    f.write('\n'.join([seq[i:i + 70] for i in range(0, len(seq), 70)]))
    f.write('\n' * 3)
    f.close()


def parse_snp_fasta(filename):
    """Reads a whole .fasta file with SNP info.
        (it occupies ORM, but can be handled)

        Returns:
            list of SNP objects

    """
    with open(filename, encoding='utf-8') as f:
        snp_fasta = f.read()
        snp_fasta = snp_fasta.split('>')[1:]        # snp_fasta[0] appears to be ''
        snp_fasta = list(sorted(set(snp_fasta)))    # remove duplicates
        snp_list = []
        for snp_fasta_item in snp_fasta:
            snp_list.append(SNP(snp_fasta_item))
    return snp_list


COMPL = {
    'A': 'T',
    'G': 'C',
    'T': 'A',
    'C': 'G',
    'N': 'N'
}


def reverse_complement(seq):
    """builds reverse complement sequence
    currently only A, T, G, C are supported"""
    return ''.join(reversed([COMPL[base] for base in seq]))


def find_snp_pos(seq, snp):
    """
    Finds SNP location based on flank overlaps.
    If the flanking sequences are from non-coding strand,
    then this function calls for reverse_complement(),
    changes flanking sequences and makes a second attempt.
    For best results, seq variable should be longer than flanks.

    :returns
      > SNP position in the sequence
        (to be used as index)
      > list of possible positions (if more than one)
    """

    seq = seq.upper()
    # print('******************************', seq) #####
    # inner variable. Letter case in original sequense remains unchanged
    # case unification allows detection of overlaps with sequences
    # that already contain different cases.
    possible_snp_pos = [0, 0, ''] # pos, probability, match diagram e.g. *** ***@**** *,
    # print('5-flank', snp.flank_5)########
    # print('3-flank', snp.flank_3)########
    # print(snp.alleles)
    for i in range(len(seq)):
        flank_5_trnc = snp.flank_5[ max(0, len(snp.flank_5)-i) : ]
        flank_3_trnc = snp.flank_3[ : min(len(seq)-i-1, len(snp.flank_3)) ]
        # 5' flank should be truncated if searching near the start
        # 3' flank should be truncated if searching near the end
        q_subseq = flank_5_trnc + "N" + flank_3_trnc
        q_subseq = q_subseq.upper() # just in case ##############
        # let it be just N between 2 joined flanks :)
        # print('88888888888888888888') #######
        # print(q_subseq) ############
        matched_nucl = 0
        match_diagram = '' # # #
        for j in range(len(q_subseq)):
            if match(seq[i - len(flank_5_trnc) + j], q_subseq[j]):
                matched_nucl += 1
                if q_subseq[j] == 'N': # 'centre' # # #
                    match_diagram += '@' # mark 'centre' # # #
                else: # # #
                    match_diagram += '*' # # #
            else: # # #
                match_diagram += ' ' # # #
        flank_match_degree = matched_nucl / len(q_subseq)
        if flank_match_degree > possible_snp_pos[-2]:  # last element is probability
            possible_snp_pos = [i, flank_match_degree, match_diagram]
        elif flank_match_degree == possible_snp_pos[-2]:
            possible_snp_pos.insert(-2, i)  # add new pos with the same probability (highly unlikely)
    # print('possible_snp_pos:', possible_snp_pos)
    if possible_snp_pos[-2] != 1:
        # print ("trying again")  ######
        # Maybe flanking sequences were from different strand.
        # Let's reverse-complement them! Maybe this will give better results
        alt_flank_5 = reverse_complement(snp.flank_3)
        alt_flank_3 = reverse_complement(snp.flank_5)
        # print('alt_5-flank', alt_flank_5)########
        # print('alt_3-flank', alt_flank_3)########
        # ...and launch the seach cycle again!
        for i in range(len(seq)):
            flank_5_trnc = alt_flank_5[ max(0, len(alt_flank_5)-i) : ]
            flank_3_trnc = alt_flank_3[ : min(len(seq)-i-1, len(alt_flank_3)) ]
            q_subseq = flank_5_trnc + "N" + flank_3_trnc
            q_subseq = q_subseq.upper()  # just in case ##############
            matched_nucl = 0
            match_diagram = ''  # # #
            for j in range(len(q_subseq)):
                if match(seq[i - len(flank_5_trnc) + j], q_subseq[j]):
                    matched_nucl += 1
                    if q_subseq[j] == 'N':  # 'centre' # # #
                        match_diagram += '@'  # mark 'centre' # # #
                    else:  # # #
                        match_diagram += '*'  # # #
                else:  # # #
                    match_diagram += '.'  # # #
            flank_match_degree = matched_nucl / len(q_subseq)
            if flank_match_degree > possible_snp_pos[-2]:  # last element is probability
                possible_snp_pos = [i, flank_match_degree, match_diagram]
                snp.flank_5 = alt_flank_5  # reverse complement flanks work better
                snp.flank_3 = alt_flank_3  # so let's reassign them
                snp.flanks_were_reversed = True  # now we know that these flanks weren't from original snp fasta file
            elif flank_match_degree == possible_snp_pos[-2]:
                possible_snp_pos.insert(-2, i)  # add new pos with the same probability (highly unlikely)
    print('possible_snp_pos:', possible_snp_pos)
    return possible_snp_pos
    # if len(possible_snp_pos) > 2 than more than one position is possible. Please check manually.


def find_snp_chrpos(rs):
    """Finds SNP chromosome position
    based on rs number
    If it's present in the cache file
    the function takes it from the file
    If not, the function makes a request

    This function works with its file autonomously
    Abscence of this file won't raise error
    Deleting this file won't raise error
    But unexpected content in the file may raise error
    Don't edit this file manually.

    :argument
        rs number with 'rs' (str)
    :returns
        chrpos (int)
    """
    with open('snpchrpos_cache.snpchrpos', 'at') as f:
        pass  # creates file if it doesn't exist
    with open('snpchrpos_cache.snpchrpos', 'rt') as f:
        for line in f:
            if line.startswith(rs):
                chrpos = int(line.strip().split('\t')[1])
                return chrpos
    # if nothing was found in dump file - let's search in web.
    URL = "https://www.ncbi.nlm.nih.gov/snp/?term=" + str(rs)
    res = requests.get(URL, timeout=300)
    if res.status_code == 200:
        try:
            chrpos = re.search(
                r"""Chromosome: </dt><dd>\d*?\:(\d*?)<br""",
                res.text.replace("\n", "")).groups()[0]
            with open('snpchrpos_cache.snpchrpos', 'at') as f:
                f.write(rs + '\t' + str(chrpos) + '\n')
        except AttributeError: # couldn't find chrpos in html file (unlikely)
            return
        return chrpos


def mark_snp(cons_cdna_seq, snp):
    """Marks given SNP in cons_cdna_seq
    using letters for ambiguous nucleotides
    :returns:
        sequence containing uppercase letter in SNP site
        indicating ambiguous nucleotide
    """
    print(snp.uid, end='\t')  ######
    if snp.clas == "snp":  # not an indel
        # print("snp.clas == snp") ######
        snp.position = find_snp_pos(cons_cdna_seq, snp)
        # print('snp.position:', snp.position)
        if len(snp.position) == 3:  # 3 elements - one best position, flank match degree, flank match diagram
            # print('cons_cdna_seq[snp.position]', cons_cdna_seq[snp.position[0]]) #########
            # print('snp.alleles', snp.alleles) #########
            # print('snp.flanks_were_reversed', snp.flanks_were_reversed) ###########
            if not snp.flanks_were_reversed:
                # take nucl from snp.alleles, find respective cons. letter and use it as a mark
                print('case1. snp alleles:', snp.alleles, 'nucl in original: ',
                      cons_cdna_seq[snp.position[0]].upper())  #########
                for key in NUCL_CODE:
                    if NUCL_CODE[key] == snp.alleles:
                        print(key, 'means', NUCL_CODE[key], 'I decided to use', key)  #####
                        cons_cdna_seq = cons_cdna_seq[: snp.position[0]] + key + cons_cdna_seq[snp.position[0] + 1:]
                        break
            elif snp.flanks_were_reversed:
                # don't take nucl from snp.alleles - they are from another strand
                # Instead, find respective compl. nucl, and respective consensus letter to mark them
                print('case2. snp alleles:', snp.alleles, 'nucl in original: ',
                      cons_cdna_seq[snp.position[0]].upper())  #########
                amb_nucl_compl_strand = set()
                for nucl in list(snp.alleles):
                    amb_nucl_compl_strand.add(COMPL[nucl])
                print('strands were reversed. It means that coding strand contains', amb_nucl_compl_strand)  ####
                for key in NUCL_CODE:
                    if NUCL_CODE[key] == amb_nucl_compl_strand:
                        print(key, 'means', NUCL_CODE[key], 'I decided to use', key)  #####
                        cons_cdna_seq = cons_cdna_seq[: snp.position[0]] + key + cons_cdna_seq[snp.position[0] + 1:]
                        break
                else:
                    pass  # can't mark SNP. Something doesn't coincide. It's better not to spoil the sequence.
    print('-------------------------------------')
    return cons_cdna_seq


# def is_self_sufficient(original_seq, snp):
#     """
#     This function can check whether a particular SNP
#     can cause ARE formation without any other SNPs
#     Works only for well located SNPs (not with 2+ possible positions)
#     :returns: tuple:
#       > is self sufficient? (True or False) (boolean)
#       > message (str)
#     """
#
#     original_seq = original_seq.lower()
#     cons_seq = mark_snp(original_seq, snp)
#     alleles = []
#     for nucleotide in NUCL_CODE[cons_seq[snp.position[0]]] - set('U'):
#         # for all the nucleotides in SNP position
#         allele = snp.uid + '_' + nucleotide
#         full_seq = original_seq[: snp.position[0]] + nucleotide + original_seq[snp.position[0] + 1:].lower()
#         print(type(full_seq))
#         print('000000000', full_seq)
#         for ARE in AREs:
#             full_seq = scan_for_overlaps(full_seq, ARE, 0)
#         region_seq = full_seq[max(0, snp.position[0] - 25) : min(snp.position[0] + 26, len(full_seq))]
#         if full_seq[snp.position[0]].upper() == original_seq[snp.position[0]].upper():
#             allele += ' (original sequence)'
#         contains_are = full_seq[snp.position[0]].isupper()
#         if contains_are:
#             some_alleles_contain_are = True
#         if ' (original sequence)' in allele and not contains_are:
#             original_seq_without_are = True
#         alleles.append({
#             'allele': allele,
#             'nucl': nucleotide,
#             'region_seq': region_seq,
#             'ARE': contains_are,
#         })
#     is_sufficient = original_seq_without_are and some_alleles_contain_are
#     return is_sufficient, alleles


def is_self_sufficient(original_seq, snp):
    """
    This function can check whether a particular SNP
    can cause ARE formation without any other SNPs
    Works only for well located SNPs (not with 2+ possible positions)
    Doesn't work for indels
    :returns: tuple:
      > is self sufficient? (True or False) (boolean)
      > message (str)
    """
    if snp.clas == 'in-del':
        return False, "it's in-del, can't check it"
    original_seq = original_seq.lower()
    cons_seq = mark_snp(original_seq, snp).lower()
    for ARE in AREs:
        cons_seq = scan_for_overlaps(cons_seq, ARE)[0]
        original_seq = scan_for_overlaps(original_seq, ARE)[0]
    if cons_seq[snp.position[0]].isupper() and original_seq[snp.position[0]].islower():
        title = 'SNP ' + snp.uid
        or_seq_fragm = original_seq[max(0, snp.position[0] - 25) : min(snp.position[0] + 26, len(original_seq))]
        cons_seq_fragm = cons_seq[max(0, snp.position[0] - 25) : min(snp.position[0] + 26, len(cons_seq))]
        footer = cons_seq[snp.position[0]].upper() + ' means ' + ' or '.join(list(NUCL_CODE[cons_seq[snp.position[0]].upper()]))
        message = '# ' + '\n# '.join([title, or_seq_fragm, cons_seq_fragm, footer])
        return True, message
    else:
        return False, ''


def _main():

    # 1. Input from .gb file
    ### filename1 = input('Enter the filename (without extension)\n')
    filename1 = "hKLF10"
    DEF, ACC, ORG, prot_seq, cdna_seq = read_from_GenBank_file(filename1)
    cons_cod_seq = reverse_translate(prot_seq)

    # 2. Processing cdna_seq
    cdna_seq, _3UTRstart = scan_for_overlaps(cdna_seq, cons_cod_seq)
    for ARE in AREs:
        cdna_seq = scan_for_overlaps(cdna_seq, ARE, _3UTRstart)[0]

    # 3. Output of cdna_seq with highlighted AREs and
    print('''Enter the name of the output file for original cdna sequence(without extension).
            * If the file exists, the sequence will be appended to it.
            * If the file doesn't exist, it will be created.''')
    filename2 = input()
    # print(cdna_seq)
    output_as_FASTA(filename2, cdna_seq, ACC, ORG)  # species name only

    # 4. Input SNP info from SNP .fasta file. Generate list of SNP objects
    ### filename3 = input('Enter the filename (with extension)\n')
    filename3 = 'snp_result.txt'
    snp_list = parse_snp_fasta(filename3)

    # 5. Obtaining chrpos for all SNPs
    for snp in snp_list:
        snp.chrpos = find_snp_chrpos(snp.uid)
        print(snp.uid, '\tchrpos:', snp.chrpos) ######


    # 6. Mark SNPs in cons_cdna_seq
    cons_cdna_seq = cdna_seq
    for snp in snp_list:
        cons_cdna_seq = mark_snp(cons_cdna_seq, snp)
    cons_cdna_seq = cons_cdna_seq.lower()

    # 7. Compare snp chrpos with assumed position in this sequence
    for snp in snp_list:
        try:
            print(snp.uid, snp.position[0], snp.chrpos, snp.chrpos+snp.position[0])
        except TypeError:
            print(snp.uid, snp.clas) # no position determined
    # snp.chrpos + snp.position[0] is always the same
    # this means that all snp positions determined by find_snp_pos were correct (at least for hKLF10 3'UTR)
    # this script uses 2 different methods to predict and validate SPN position in the sequence
    # it' a very powerful tool. In order to use 100% of its abilities it's recommended to
    # to improve comparison algorithm.
    # thois should be done before marking these SNPs - i.e. mark_snps function should call validation function

    # 8. mark CDS and AREs in consensus sequence and output it side-by-side with original sequence
    # 8.1. Processing cons_cdna_seq
    cons_cdna_seq, _3UTRstart = scan_for_overlaps(cons_cdna_seq, cons_cod_seq)
    for ARE in AREs:
        cons_cdna_seq = scan_for_overlaps(cons_cdna_seq, ARE, _3UTRstart)[0]

    # 8.2. Output of cdna_seq with highlighted AREs and
    print('''Enter the name of the output file for consensus cna sequence (without extension).
            * If the file exists, the sequence will be appended to it.
            * If the file doesn't exist, it will be created.''')
    filename2 = input()
    output_as_FASTA(filename2, cons_cdna_seq, ACC, ORG)  # species name only

    print('Original cDNA sequence vs consensus cDNA sequence')
    print(cdna_seq)
    print(cons_cdna_seq)

    # 9. Now we have consensus sequence. We can see that some new AREs appeared
    # But many of these AREs contain several SNPs
    # Some of them, though, contain only one SNP
    # It would be useful to find SNPs that can cause ARE without any other SNPs
    for snp in snp_list:
        res = is_self_sufficient(cdna_seq.lower(), snp)
        if res[0]:
            print(res[1])


if __name__ == "__main__":
    _main()