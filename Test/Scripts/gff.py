# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 14:41:15 2024

@author: Anna
"""

#extract the loci which are genes/cds/etc from the gff3 file
# =============================================================================
# def GFFExtractLoci(gff_input, searchword):
#     import re
#     gff_feature = list()
#     i = 0
#     while i < len(gff_input):
#         x = re.findall(searchword, gff_input[i])
#         if x:
#             gff_feature.append(gff_input[i].split("\t"))
#         i=i+1
#     return gff_feature
#         
# =============================================================================

def GFFExtractLoci(gff_input, searchword):
    gff_feature = []
    append = gff_feature.append
    sw = searchword  # lokale Variable
    
    for line in gff_input:
        if line.find(sw) != -1:
            append(line.split("\t"))

    return gff_feature
