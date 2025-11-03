# -*- coding: utf-8 -*-
#Reformat the fasta
def FastaCleanup(fasta_list):
    import re #import the regex library
    i=0
    fa = ["",""] # An empty list, which is needed for the fasta 
    while i < len(fasta_list):
        x = re.findall("\A>", fasta_list[i])
        if x:
            fa[0]=fasta_list[i]
        else:
            fa[1]=fa[1]+fasta_list[i]
        i=i+1
    fa[1]=''.join(fa[1].split('\n'))
    return fa

def ExtractRegion(lower, upper, direction, data):
    import re
    subregion = str()
    for i in data:
        if direction == '-':
            if not re.search("^>", i):
                subregion = i[lower-1:upper-1]
        elif direction == '+':
            if not re.search("^>", i):
                subregion_tmp = i[lower-1:upper-1]
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
