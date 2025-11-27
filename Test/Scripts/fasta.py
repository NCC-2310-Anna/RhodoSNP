# -*- coding: utf-8 -*-
#Reformat the fasta
def FastaCleanup(fasta_list):
    headers = []
    seq_chunks = []  # statt seqs = [] mit +=

    h_append = headers.append
    c_append = seq_chunks.append

    for line in fasta_list:
        line = line.strip()
        if line.startswith(">"):
            h_append(line)
            c_append([])  # neue Liste für seq lines
        else:
            if not seq_chunks:
                h_append(None)
                c_append([])
            seq_chunks[-1].append(line)

    # Jetzt joinen wir die Sequenzen einmal
    seqs = ["".join(ch) for ch in seq_chunks]
    return list(zip(headers, seqs))

def ExtractRegion(lower, upper, direction, data):
    subregion = str()
    if direction == '-':
        subregion = data[lower-1:upper-1]
    elif direction == '+':
        subregion_tmp = data[lower-1:upper-1]
        subregion = list()
        for j in subregion_tmp:
            subregion.append(j)
        subregion = ''.join(subregion)
    return subregion


def ReverseString(data):
    data = data[::-1]
    return data

def KmerSplit(data, kmer_size):
    split_string_list = [data[x:x+kmer_size] for x in range(0,len(data),kmer_size)]
    return split_string_list
