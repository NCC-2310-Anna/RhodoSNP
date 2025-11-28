# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 14:06:51 2024

@author: Anna
"""
def merge_vcf(vcf1, vcf2, vcf3, searchword):
    from collections import Counter

    # Build lookup dictionaries for faster access
    vcf_dicts = [{rec[0]: rec for rec in vcf} for vcf in (vcf1, vcf2, vcf3)]

    # Count all positions across files
    all_positions = sum((list(d.keys()) for d in vcf_dicts), [])
    pos_counts = Counter(all_positions)

    # Keep only positions present in at least 2 VCFs
    common_positions = [p for p, c in pos_counts.items() if c >= 2]

    def fill_missing(i, j, k):
        """If INFO are '.', take from next available VCF (priority i → j → k)"""
        info = i[3] if len(i) > 3 and i[3] != '.' else (j[3] if len(j) > 3 and j[3] != '.' else (k[3] if len(k) > 3 else '.'))
        return info

    merged_records = []
    for pos in common_positions:
        a = vcf_dicts[0].get(pos, [pos, '.', '.', '.'])
        b = vcf_dicts[1].get(pos, [pos, '.', '.', '.'])
        c = vcf_dicts[2].get(pos, [pos, '.', '.', '.'])

        # Fill missing fields by borrowing from another VCF
        i_info = fill_missing(a, b, c)

        # Store merged line
  #      merged_records.append([pos, i_ref, i_alt, j_ref, j_alt, k_ref, k_alt, i_info])
        merged_records.append(
            [pos, a[1], a[2], b[1], b[2], c[1], c[2], i_info]
        )

    # Robust SNP/Indel detection
    def is_snp(record):
        ref1, alt1, ref2, alt2, ref3, alt3 = record[1:7]

        # Flatten ALT alleles if comma separated
        alts = []
        for allele in (alt1, alt2, alt3):
            if allele == '.':
                alts.append('.')
            else:
                alts.extend(allele.split(','))

        return all(len(r) == 1 or r == '.' for r in (ref1, ref2, ref3)) and \
               all(len(a) == 1 or a == '.' for a in alts)

    # Return filtered results
    searchword = searchword.lower()
    if searchword == "all":
        return merged_records
    if searchword == "snp":
        return [r for r in merged_records if is_snp(r)]
    if searchword == "indel":
        return [r for r in merged_records if not is_snp(r)]

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