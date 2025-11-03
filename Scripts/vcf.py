# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 14:06:51 2024

@author: Anna
"""

def MergeVCF(vcf1, vcf2, vcf3, searchword):
    from collections import Counter

    # Use Dictionaries for faster Data accesss
    vcf_dicts = [{rec[0]: rec for rec in vcf} for vcf in [vcf1, vcf2, vcf3]]
    
    # collect all positions
    positions = []
    for d in vcf_dicts:
        positions.extend(d.keys())
    
    # Positionen zählen
    pos_counts = Counter(positions)
    
    # Nur Positionen mit mindestens 2 Vorkommen
    common_positions = [pos for pos, count in pos_counts.items() if count >= 2]
    
    merged_records = []
    for pos in common_positions:
        # Hole Werte, wenn vorhanden, sonst Platzhalter
        i = vcf_dicts[0].get(pos, [pos, '.', '.', '.', '.'])
        j = vcf_dicts[1].get(pos, [pos, '.', '.', '.', '.'])
        k = vcf_dicts[2].get(pos, [pos, '.', '.', '.', '.'])
        
        # Format: Pos, i_REF, i_ALT, j_REF, j_ALT, k_REF, k_ALT, INFO
        merged = [pos, i[1], i[2], j[1], j[2], k[1], k[2], i[3] if len(i) > 3 else '.']
        merged_records.append(merged)

    def is_snp(record):
        alleles = record[1:7]
        # ein SNP = Länge aller Einträge == 1 oder '.'
        return all(len(a) == 1 or a == '.' for a in alleles)

    if searchword == "all":
        return merged_records
    elif searchword == "indel":
        return [rec for rec in merged_records if not is_snp(rec)]
    elif searchword == "snp":
        return [rec for rec in merged_records if is_snp(rec)]
    else:
        print("Searchword must be 'snp', 'indel' oder 'all'!")
        return []


#extract a list of SNPs from the vcf File 
def VCFExtractSNP(vcf_input):
    import re
    result_1 = []
    result_2 = []
    Chrom = []
    SNP = []
    Ref = []
    Alt = []
    SNP_Ind = 0
    Alt_Ind = 0
    Ref_Ind = 0
    Chr_Ind = 0
    for x in vcf_input:
        i = re.search("(\A##{1})",x)
        a = re.search("(\A#{1})", x)
        if not i:
            result_1.append(x)
        if a:
            c =0
            while c < len(x.split("\t")):
                if x.split("\t")[c] == "CHROM":
                    Chr_Ind = c
                elif x.split("\t")[c] == "POS":
                    SNP_Ind = c
                elif x.split("\t")[c] == "REF":
                    Ref_Ind = c
                elif x.split("\t")[c] == "ALT":
                    Alt_Ind = c
                c = c+1
    for x in result_1:
        Chrom.append(x.split('\t')[Chr_Ind])
        SNP.append(x.split('\t')[SNP_Ind])
        Ref.append(x.split('\t')[Ref_Ind])
        Alt.append(x.split('\t')[Alt_Ind])
    i =0
    while i < len(SNP):
        tmp_list = [SNP[i],Ref[i],Alt[i], Chrom[i]]
        result_2.append(tmp_list)
        i=i+1
    return result_2

def SplitVCF(data,searchword):
    kind_mutation = list()
    if searchword == "indel":
        for i in data:
            if (len(i[1])!=1 or len(i[2])!=1):
                kind_mutation.append(i)
    elif searchword == "snp":
        kind_mutation.append(data[0])
        for i in data:
            if (len(i[1])==1 and len(i[2])==1):
                kind_mutation.append(i)
    else:
        print("Searchword must be 'snp' or 'indel'!")
    return kind_mutation

#extract all SNPs from the truth dataset provided by bush et al
def TruthExtractSNP(tsv_input):
    import re
    result1 = []
    for x in tsv_input:
        a = re.search("\AChr", x)
        if not a:
            result1.append(x.split('\t'))
    return result1

def split_vcf_into_blocks(input_vcf, block_size=10000, pos_index = 0):
    blocks = {}
    current_headers = []
    for line in input_vcf:
        if line[pos_index] != "POS":
            pos = int(line[pos_index])
            block_start = (pos // block_size) * block_size + 1
            block_end = block_start + block_size - 1
            block_key = f"{block_start}-{block_end}"

            if block_key not in blocks:
                blocks[block_key] = current_headers.copy()  # Header pro Block

            blocks[block_key].append(line)

    return blocks