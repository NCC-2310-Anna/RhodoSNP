# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 14:06:51 2024

@author: Anna
"""

def MergeVCF(vcf1, vcf2, vcf3, searchword):
    vcf_all =[] 
    vcf_common = []
    counter = 0
    for i in vcf1:
        for j in vcf2:
            for k in vcf3:
                if i[0] == j[0] and i[0] == k[0]:
                    vcf_all.append([i[0],i[1],i[2],j[1],j[2],k[1],k[2],i[3]])
                    counter = counter + 1
    if searchword == "all":
        return vcf_all
    elif searchword == "indel":
        for i in vcf_all:
            if (len(i[1]) != 1 or len(i[2]) != 1 or\
                len(i[3]) != 1 or len(i[4]) != 1 or\
                len(i[5]) != 1 or len(i[6]) != 1):
                vcf_common.append(i)
        vcf_all = vcf_common
        return vcf_all
    
    elif searchword == "snp":
        vcf_common.append(vcf_all[0])
        for i in vcf_all:
            if (len(i[1])==1 and len(i[2])==1 and\
                len(i[3])==1 and len(i[4])==1 and\
                len(i[5])==1 and len(i[6])==1):
                vcf_common.append(i)
    else:
        print("Searchword must be 'snp', 'indel' or 'all'!")
    return vcf_common

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