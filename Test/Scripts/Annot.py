#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Variant-to-Protein Annotation Script
Optimized version for pipeline integration
Author: annan
Date: 2025-10-27
"""

import sys, time
from argparse import ArgumentParser

sys.path.append("PythonScripts")
import ReadFiles as Rf
import vcf, fasta, gff, translation_table


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main():
    parser = ArgumentParser(description="Merge VCFs and generate protein-level SNP annotation.")
    parser.add_argument("-f", "--fasta", required=True, help="Path to reference FASTA file")
    parser.add_argument("-g", "--gff", required=True, help="Path to GFF genome annotation file")
    parser.add_argument("-v1", "--vcf1", required=True, help="VCF from bcftools")
    parser.add_argument("-v2", "--vcf2", required=True, help="VCF from freebayes")
    parser.add_argument("-v3", "--vcf3", required=True, help="VCF from lofreq")
    parser.add_argument("-o1", "--outputProtein", required=True, help="Output: mutated protein FASTA")
    parser.add_argument("-o2", "--outputAnnotation", required=True, help="Output: amino acid change list")
    parser.add_argument("-o3", "--outputSNP", required=True, help="Output: consensus SNP list")
    args = parser.parse_args()
    
        # --- Load Input Files ---
    # log("Loading input files ...")
    # fasta_file = Rf.ReadFile(args.fasta)
    # gff_file = Rf.ReadFile(args.gff)
    # vcf_bcf = Rf.ReadFile(args.vcf1)
    # vcf_fb = Rf.ReadFile(args.vcf2)
    # vcf_lf = Rf.ReadFile(args.vcf3)
    
        # --- Load Input Files ---
    log("Loading input files ...")
    fasta_file = Rf.ReadFile("/home/anna-nauruschat/Dokumente/Phd/Projects/Targets/Cereibacter_sphaeroides_reoriented.fasta")
    gff_file = Rf.ReadFile("/home/anna-nauruschat/Dokumente/Phd/Projects/Targets/annot.gff")
    vcf_bcf = Rf.ReadFile("/home/anna-nauruschat/Dokumente/Phd/Projects/Targets/SRR5769006BC.vcf")
    vcf_fb = Rf.ReadFile("/home/anna-nauruschat/Dokumente/Phd/Projects/Targets/SRR5769006FB.vcf")
    vcf_lf = Rf.ReadFile("/home/anna-nauruschat/Dokumente/Phd/Projects/Targets/SRR5769006LF.vcf")
    
    # --- Prepare data ---
    log("Cleaning FASTA and extracting coding regions ...")
    fasta_new = fasta.FastaCleanup(fasta_file)
    gff_loci = gff.GFFExtractLoci(gff_file, "CDS")
    
    log("Extracting SNPs ...")
    vcf_new_bcf = vcf.VCFExtractSNP(vcf_bcf)
    vcf_new_fb = vcf.VCFExtractSNP(vcf_fb)
    vcf_new_lf = vcf.VCFExtractSNP(vcf_lf)
    
    # --- Merge VCFs ---
    log("Merging SNPs from all callers ...")
    if len(vcf_new_bcf) < 1000 or len(vcf_new_fb) < 1000 or len(vcf_new_lf) < 1000:
        vcf_snp = vcf.merge_vcf(vcf_new_bcf, vcf_new_fb, vcf_new_lf, "all")
    else:
        log("Large files detected — splitting and merging by blocks ...")
        blocks = [
            vcf.split_vcf_into_blocks(vcf_new_bcf),
            vcf.split_vcf_into_blocks(vcf_new_fb),
            vcf.split_vcf_into_blocks(vcf_new_lf),
        ]
        vcf_snp = []
        for key in set(blocks[0].keys()) & set(blocks[1].keys()) & set(blocks[2].keys()):
            vcf_snp.extend(vcf.merge_vcf(blocks[0][key], blocks[1][key], blocks[2][key], "all"))
    log(f"Consensus SNPs found: {len(vcf_snp)-1}")
    
    # --- Consensus SNPs cleanup ---
    vcf_snp_2 = vcf.ConflictDetection(vcf_snp)
    vcf_snp = vcf_snp_2
    
    # Remove the variables we dont need anymore...
    del vcf_bcf, vcf_fb, vcf_lf, vcf_snp_2, vcf_new_bcf, vcf_new_fb, vcf_new_lf, fasta_file, gff_file
    
    # --- Extract coding regions ---
    log("Extracting CDS subregions ...")
    
    # Header row
    subregions = [["Lower Limit", "Upper Limit", "Direction", "Info", "Region", "Chromosome"]]
    
    # lokale Bindungen für mehr Speed
    ex_region = fasta.ExtractRegion
    fa_new = fasta_new
    sr_append = subregions.append
    
    for i in gff_loci:
        lower = int(i[3])
        upper = int(i[4]) + 1  # nur einmal berechnen
        direction = i[6]
        info = i[8]
        chrom = int(i[0])
    
        seq = fa_new[int(chrom)-1][1]  # Chromosom lookup nur 1× pro Loop
        region = ex_region(lower, upper, direction, seq)
    
        sr_append([lower, upper, direction, info, region, chrom])
        
    
    # --- Map SNPs to subregions ---
    log("Mapping SNPs to CDS regions ...")
    subregion_with_snp = []
    for region in subregions[1:]:
        for snp in vcf_snp:
            if snp[0].isdigit():
                if region[0] <= int(snp[0]) <= region[1] and region[5]==int(snp[3]):
                    subregion_with_snp.append(
                        [snp[0], snp[1], snp[2], *region[:4], region[4], snp[3]]
                    )
    
    # --- Generate mutated FASTA ---
    log("Generating mutated DNA sequences ...")
    multifasta = []
    for entry in subregion_with_snp:
        pos, ref, alt, start, end, direction, info, seq, chrom = entry
        if(len(ref)==1 and len(alt)==1):
            seq = list(seq)
            idx = int(pos) - int(start)
            if 0 <= idx < len(seq) and seq[idx] == ref:
                seq[idx] = alt
                mutated = "".join(seq)
                if direction == "-":
                    mutated = fasta.ReverseString(mutated).translate(str.maketrans("ATCG", "TAGC"))
                    multifasta.append((f">{pos}_{ref}_{alt}_{start}_{end}_{direction}", mutated))
            else:
                log(f"[WARN] Mismatch at SNP {pos}: expected {ref}, found {alt}")
    
    # --- Translate DNA to protein ---
    log("Translating DNA to protein ...")
    protein_records = []
    for header, seq in multifasta:
        prot = translation_table.translateDNA(fasta.KmerSplit(seq, 3))
        protein_records.append((header, prot))
        
    # --- Build annotation table ---
    log("Building amino acid annotation ...")
    annotation = []
    for header, prot_seq in protein_records:
        parts = header.split("_")
        pos, ref, alt, start, end, direction = parts[0:6]
        for region in subregions[1:]:
            if int(start) == region[0]:
                info = region[3].replace("\n", "").split(";")
                annotation.append(
                    "\t".join([
                        pos,
                        info[-3] if len(info) >= 3 else "",
                        info[-2] if len(info) >= 2 else "",
                        ref, alt, direction
                    ])
                )
    
    # --- Write output files ---
    log("Writing output files ...")
    with open(args.outputAnnotation, "w") as f:
        f.write("Pos\tProduct\tProtein\tRef\tAlt\tStrand\n")
        f.write("\n".join(annotation))
    
    with open(args.outputProtein, "w") as f:
        for h, s in protein_records:
            f.write(f"{h}\n{s}\n")
    
    
    with open(args.outputSNP, "w") as f:
        for snp in vcf_snp:
            f.write("\t".join(map(str, snp)) + "\n")
    
    log("All done!")


if __name__ == "__main__":
    main()