from pathlib import Path

def analysis_of_DNA(name):
    aminoAcids = {
        'Met': ['AUG'],
        'Phe': ['UUU', 'UUC'],
        'Leu': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
        'Ile': ['AUU', 'AUC', 'AUA'],
        'Val': ['GUU', 'GUC', 'GUA', 'GUG'],
        'Ser': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
        'Pro': ['CCU', 'CCC', 'CCA', 'CCG'],
        'Thr': ['ACU', 'ACC', 'ACA', 'ACG'],
        'Ala': ['GCU', 'GCC', 'GCA', 'GCG'],
        'Tyr': ['UAU', 'UAC'],
        'His': ['CAU', 'CAC'],
        'Gln': ['CAA', 'CAG'],
        'Asn': ['AAU', 'AAC'],
        'Lys': ['AAA', 'AAG'],
        'Asp': ['GAU', 'GAC'],
        'Glu': ['GAA', 'GAG'],
        'Cys': ['UGU', 'UGC'],
        'Trp': ['UGG'],
        'Arg': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'Gly': ['GGU', 'GGC', 'GGA', 'GGG'],
        'Stop': ['UAA', 'UAG', 'UGA']
    }

    RNA_seq_str = ''

    aminoAcidSeq = []

    AA_HM_dict = {}

    AA_HM_all = 0

    with open(name, mode='r', encoding='utf-8') as DNA_seq:
        for line in DNA_seq:
            for letter in line:
                if letter == 'A':
                    RNA_seq_str += 'U'
                elif letter == 'T':
                    RNA_seq_str += 'A'
                elif letter == 'C':
                    RNA_seq_str += 'G'
                elif letter == 'G':
                    RNA_seq_str += 'C'
                else:
                    RNA_seq_str += letter

    name = name[:-4]

    with open(f'{name}_butRNA.txt', mode='w', encoding='utf-8') as RNA_seq:
        RNA_seq.write(RNA_seq_str)

    with open(f'{name}_butRNA.txt', mode='r', encoding='utf-8') as RNA_seq_trans:
        for line in RNA_seq_trans:
            for index in range(0, len(line), 3):
                codon = line[index:index+3]
                for aminoAcidName in aminoAcids:
                    for aminoAcid in range(len(aminoAcids[aminoAcidName])):
                        if codon == aminoAcids[aminoAcidName][aminoAcid]:
                            aminoAcidSeq.append(aminoAcidName)

    with open(f'{name}_butAA.txt', mode='w', encoding='utf-8') as AA_seq:
        for AA in aminoAcidSeq:
            AA_seq.write(f'{AA}\n')

    with open(f'{name}_butAA.txt', mode='r', encoding='utf-8') as AA_seq:
        for AA in AA_seq:
            AA = AA[:-1]
            if AA not in AA_HM_dict:
                AA_HM_dict[AA] = 1
                AA_HM_all += 1
            else:
                AA_HM_dict[AA] += 1
                AA_HM_all += 1

    with open(f'{name}_butAA_butHM.txt', mode='w', encoding='utf-8') as AA_HM:
        for AAname in AA_HM_dict:
            AA_HM.write(f'{AAname}: {AA_HM_dict[AAname]} | {round((AA_HM_dict[AAname]/AA_HM_all)*100, 2)}%\n')

analysis_of_DNA('european_rabbit.txt')