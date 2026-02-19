# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 13:47:36 2026

@author: annan
"""

# Usage: QualityMetricsCalculation.py -t /path/to/truth/data -f /path/to/calling/results
from argparse import ArgumentParser

# read the files
def ReadFile(file):
     fobj_file = open(file, "r")
     file = fobj_file.readlines()
     fobj_file.close()
     return file

# extract all SNPs from the truth dataset provided by bush et al
def TruthExtractSNP(tsv_input):
    import re
    result = []
    for x in tsv_input:
        a = re.search("\AChr", x)
        if not a:
            result.append(x.split('\t'))
    return result

# extract all SNPs from the Caller TSV
def TSVExtractSNP(tsv_input):
    import re
    result = []
    for x in tsv_input:
        a = re.search("\APOS", x)
        if not a:
            result.append(x.split('\t'))
    return result

# convert a list of strings with a tab separator into a list of lists 
def ConvertToDict(inputData, keyIndex):
    outputDict = []
    for line in inputData:
        temp=line.split('\t')
        tempDict = {temp[keyIndex]: temp[:keyIndex]+temp[keyIndex+1:]}
        outputDict.append(tempDict)
    return outputDict

def CalculateQualityMetrics(truthData, snpData):
    truePositive = 0
    falsePositive = 0
    falseNegative = 0

    snpKeys = {k for d in snpData for k in d}
    truthKeys = {k for d in truthData for k in d}

    for snp_dict in snpData:
        for pos, values in snp_dict.items():
            if pos in truthKeys:
                # weitere Bedingungen prüfen
                truth_values = next(d[pos] for d in truthData if pos in d)
                if values[0] == truth_values[0] and values[3]=="SNP\n" and truth_values[4]=="*\n":  # z.B. REF und ALT müssen übereinstimmen
                    truePositive += 1
            else:
                if values[3]=="SNP\n":  # z.B. REF und ALT müssen übereinstimmen
                    falsePositive += 1

    for pos in truthKeys:
        if pos not in snpKeys:
            falseNegative += 1
            
    precision = (truePositive)/(truePositive+falsePositive)
    recall = (truePositive)/(truePositive+falseNegative)
    f1Score = (2*precision*recall)/(precision+recall)
    return [truePositive, falsePositive, falseNegative, precision, recall, f1Score]

# get metadata (caller) from the filename
def ParseMetadata(filePath):
    pathList = filePath.split("/")
    fileName = pathList[len(pathList)-1]
    del pathList
    fileTemp = fileName.replace(".tsv", "")
    metadataList = fileTemp.split("_")
    caller = metadataList[len(metadataList)-1]
    return [fileName, caller]
    
def main():
    parser = ArgumentParser(description="Calculate EvalData (TP, TN, FP, FN)")
    parser.add_argument("-t", "--truthData", required=True, help="Path to truth data file")
    parser.add_argument("-f", "--snpData", required=True,  help="Path to snp data file")
    args = parser.parse_args()
    
    truthData = ConvertToDict(ReadFile(args.truthData), 1)
    snpData = ConvertToDict(ReadFile(args.snpData), 0)
    metadata = ParseMetadata(args.snpData)
    qualMetrics = CalculateQualityMetrics(truthData, snpData)
    print(*metadata, *qualMetrics,  sep='\t')

            
main()
