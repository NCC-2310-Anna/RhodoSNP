# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 09:37:31 2024

@author: Anna
"""
# 1. Standard Code (transl_table=1)
table1 = {
    "TTT": ["F", "Standard"], "TTC": ["F", "Standard"], "TTA": ["L", "Standard"], "TTG": ["L", "Standard"],
    "CTT": ["L", "Standard"], "CTC": ["L", "Standard"], "CTA": ["L", "Standard"], "CTG": ["L", "Standard"],
    "ATT": ["I", "Standard"], "ATC": ["I", "Standard"], "ATA": ["I", "Standard"], "ATG": ["M", "start"],  # Start codon
    "GTT": ["V", "Standard"], "GTC": ["V", "Standard"], "GTA": ["V", "Standard"], "GTG": ["V", "Standard"],
    
    "TCT": ["S", "Standard"], "TCC": ["S", "Standard"], "TCA": ["S", "Standard"], "TCG": ["S", "Standard"],
    "CCT": ["P", "Standard"], "CCC": ["P", "Standard"], "CCA": ["P", "Standard"], "CCG": ["P", "Standard"],
    "ACT": ["T", "Standard"], "ACC": ["T", "Standard"], "ACA": ["T", "Standard"], "ACG": ["T", "Standard"],
    "GCT": ["A", "Standard"], "GCC": ["A", "Standard"], "GCA": ["A", "Standard"], "GCG": ["A", "Standard"],
    
    "TAT": ["Y", "Standard"], "TAC": ["Y", "Standard"], "TAA": ["*", "stop"], "TAG": ["*", "stop"], "TGA": ["*", "stop"],
    "CAT": ["H", "Standard"], "CAC": ["H", "Standard"],
    "CAA": ["Q", "Standard"], "CAG": ["Q", "Standard"],
    "AAT": ["N", "Standard"], "AAC": ["N", "Standard"],
    "AAA": ["K", "Standard"], "AAG": ["K", "Standard"],
    "GAT": ["D", "Standard"], "GAC": ["D", "Standard"],
    "GAA": ["E", "Standard"], "GAG": ["E", "Standard"],
    
    "TGT": ["C", "Standard"], "TGC": ["C", "Standard"],
    "TGG": ["W", "Standard"],
    "CGT": ["R", "Standard"], "CGC": ["R", "Standard"], "CGA": ["R", "Standard"], "CGG": ["R", "Standard"],
    "AGT": ["S", "Standard"], "AGC": ["S", "Standard"],
    "AGA": ["R", "Standard"], "AGG": ["R", "Standard"],
    "GGT": ["G", "Standard"], "GGC": ["G", "Standard"], "GGA": ["G", "Standard"], "GGG": ["G", "Standard"]
}

# 2. Vertebrate Mitochondrial Code (transl_table=2)
table2 = {
    "TTT": ["F", "Standard"], "TTC": ["F", "Standard"],
    "TTA": ["L", "Standard"], "TTG": ["L", "Standard"],
    "CTT": ["L", "Standard"], "CTC": ["L", "Standard"], "CTA": ["L", "Standard"], "CTG": ["L", "Standard"],
    "ATT": ["I", "start"], "ATC": ["I", "start"], "ATA": ["M", "start"],
    "ATG": ["M", "start"],
    "GTT": ["V", "Standard"], "GTC": ["V", "Standard"], "GTA": ["V", "Standard"], "GTG": ["V", "start"],
    
    "TCT": ["S", "Standard"], "TCC": ["S", "Standard"], "TCA": ["S", "Standard"], "TCG": ["S", "Standard"],
    "CCT": ["P", "Standard"], "CCC": ["P", "Standard"], "CCA": ["P", "Standard"], "CCG": ["P", "Standard"],
    "ACT": ["T", "Standard"], "ACC": ["T", "Standard"], "ACA": ["T", "Standard"], "ACG": ["T", "Standard"],
    "GCT": ["A", "Standard"], "GCC": ["A", "Standard"], "GCA": ["A", "Standard"], "GCG": ["A", "Standard"],
    
    "TAT": ["Y", "Standard"], "TAC": ["Y", "Standard"],
    "TAA": ["*", "stop"], "TAG": ["*", "stop"], "TGA": ["W", "Standard"],
    "CAT": ["H", "Standard"], "CAC": ["H", "Standard"],
    "CAA": ["Q", "Standard"], "CAG": ["Q", "Standard"],
    "AAT": ["N", "Standard"], "AAC": ["N", "Standard"],
    "AAA": ["K", "Standard"], "AAG": ["K", "Standard"],
    "GAT": ["D", "Standard"], "GAC": ["D", "Standard"],
    "GAA": ["E", "Standard"], "GAG": ["E", "Standard"],
    
    "TGT": ["C", "Standard"], "TGC": ["C", "Standard"],
    "TGG": ["W", "Standard"],
    "CGT": ["R", "Standard"], "CGC": ["R", "Standard"], "CGA": ["R", "Standard"], "CGG": ["R", "Standard"],
    "AGT": ["S", "Standard"], "AGC": ["S", "Standard"],
    "AGA": ["*", "stop"], "AGG": ["*", "stop"],
    "GGT": ["G", "Standard"], "GGC": ["G", "Standard"], "GGA": ["G", "Standard"], "GGG": ["G", "Standard"]
}

# 3. Yeast Mitochondrial Code (transl_table=3)
table3 = {
    "TTT": ["F", "Standard"], "TTC": ["F", "Standard"],
    "TTA": ["L", "Standard"], "TTG": ["L", "Standard"],
    "CTT": ["T", "Standard"], "CTC": ["T", "Standard"], "CTA": ["T", "Standard"], "CTG": ["T", "Standard"],
    "ATT": ["I", "start"], "ATC": ["I", "Standard"], "ATA": ["M", "start"],
    "ATG": ["M", "start"],
    "GTT": ["V", "Standard"], "GTC": ["V", "Standard"], "GTA": ["V", "Standard"], "GTG": ["V", "Standard"],
    
    "TCT": ["S", "Standard"], "TCC": ["S", "Standard"], "TCA": ["S", "Standard"], "TCG": ["S", "Standard"],
    "CCT": ["P", "Standard"], "CCC": ["P", "Standard"], "CCA": ["P", "Standard"], "CCG": ["P", "Standard"],
    "ACT": ["T", "Standard"], "ACC": ["T", "Standard"], "ACA": ["T", "Standard"], "ACG": ["T", "Standard"],
    "GCT": ["A", "Standard"], "GCC": ["A", "Standard"], "GCA": ["A", "Standard"], "GCG": ["A", "Standard"],
    
    "TAT": ["Y", "Standard"], "TAC": ["Y", "Standard"],
    "TAA": ["*", "stop"], "TAG": ["*", "stop"], "TGA": ["W", "Standard"],
    "CAT": ["H", "Standard"], "CAC": ["H", "Standard"],
    "CAA": ["Q", "Standard"], "CAG": ["Q", "Standard"],
    "AAT": ["N", "Standard"], "AAC": ["N", "Standard"],
    "AAA": ["K", "Standard"], "AAG": ["K", "Standard"],
    "GAT": ["D", "Standard"], "GAC": ["D", "Standard"],
    "GAA": ["E", "Standard"], "GAG": ["E", "Standard"],
    
    "TGT": ["C", "Standard"], "TGC": ["C", "Standard"],
    "TGG": ["W", "Standard"],
    "CGT": ["R", "Standard"], "CGC": ["R", "Standard"], "CGA": ["R", "Standard"], "CGG": ["R", "Standard"],
    "AGT": ["S", "Standard"], "AGC": ["S", "Standard"],
    "AGA": ["R", "Standard"], "AGG": ["R", "Standard"],
    "GGT": ["G", "Standard"], "GGC": ["G", "Standard"], "GGA": ["G", "Standard"], "GGG":["G","Standard"]}
# 4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (transl_table=4)
table4 = {
    "TTT": ["F", "Standard"], "TTC": ["F", "Standard"],
    "TTA": ["L", "start"], "TTG": ["L", "start"],
    "CTT": ["L", "Standard"], "CTC": ["L", "Standard"], "CTA": ["L", "Standard"], "CTG": ["L", "start"],
    "ATT": ["I", "start"], "ATC": ["I", "start"], "ATA": ["I", "start"],
    "ATG": ["M", "start"],
    "GTT": ["V", "Standard"], "GTC": ["V", "Standard"], "GTA": ["V", "Standard"], "GTG": ["V", "Standard"],
    
    "TCT": ["S", "Standard"], "TCC": ["S", "Standard"], "TCA": ["S", "Standard"], "TCG": ["S", "Standard"],
    "CCT": ["P", "Standard"], "CCC": ["P", "Standard"], "CCA": ["P", "Standard"], "CCG": ["P", "Standard"],
    "ACT": ["T", "Standard"], "ACC": ["T", "Standard"], "ACA": ["T", "Standard"], "ACG": ["T", "Standard"],
    "GCT": ["A", "Standard"], "GCC": ["A", "Standard"], "GCA": ["A", "Standard"], "GCG": ["A", "Standard"],
    
    "TAT": ["Y", "Standard"], "TAC": ["Y", "Standard"],
    "TAA": ["*", "stop"], "TAG": ["*", "stop"], "TGA": ["W", "Standard"],
    "CAT": ["H", "Standard"], "CAC": ["H", "Standard"],
    "CAA": ["Q", "Standard"], "CAG": ["Q", "Standard"],
    "AAT": ["N", "Standard"], "AAC": ["N", "Standard"],
    "AAA": ["K", "Standard"], "AAG": ["K", "Standard"],
    "GAT": ["D", "Standard"], "GAC": ["D", "Standard"],
    "GAA": ["E", "Standard"], "GAG": ["E", "Standard"],
    
    "TGT": ["C", "Standard"], "TGC": ["C", "Standard"],
    "TGG": ["W", "Standard"],
    "CGT": ["R", "Standard"], "CGC": ["R", "Standard"], "CGA": ["R", "Standard"], "CGG": ["R", "Standard"],
    "AGT": ["S", "Standard"], "AGC": ["S", "Standard"],
    "AGA": ["R", "Standard"], "AGG": ["R", "Standard"],
    "GGT": ["G", "Standard"], "GGC": ["G", "Standard"], "GGA": ["G", "Standard"], "GGG": ["G", "Standard"]
}

# 5. The Invertebrate Mitochondrial Code (transl_table=5)
table5 = {
    "TTT": ["F", "Standard"], "TTC": ["F", "Standard"],
    "TTA": ["L", "Standard"], "TTG": ["L", "start"],
    "CTT": ["L", "Standard"], "CTC": ["L", "Standard"], "CTA": ["L", "Standard"], "CTG": ["L", "Standard"],
    "ATT": ["I", "start"], "ATC": ["I", "start"], "ATA": ["M", "start"],
    "ATG": ["M", "start"],
    "GTT": ["V", "Standard"], "GTC": ["V", "Standard"], "GTA": ["V", "Standard"], "GTG": ["V", "Standard"],
    
    "TCT": ["S", "Standard"], "TCC": ["S", "Standard"], "TCA": ["S", "Standard"], "TCG": ["S", "Standard"],
    "CCT": ["P", "Standard"], "CCC": ["P", "Standard"], "CCA": ["P", "Standard"], "CCG": ["P", "Standard"],
    "ACT": ["T", "Standard"], "ACC": ["T", "Standard"], "ACA": ["T", "Standard"], "ACG": ["T", "Standard"],
    "GCT": ["A", "Standard"], "GCC": ["A", "Standard"], "GCA": ["A", "Standard"], "GCG": ["A", "Standard"],
    
    "TAT": ["Y", "Standard"], "TAC": ["Y", "Standard"],
    "TAA": ["*", "stop"], "TAG": ["*", "stop"], "TGA": ["W", "Standard"],
    "CAT": ["H", "Standard"], "CAC": ["H", "Standard"],
    "CAA": ["Q", "Standard"], "CAG": ["Q", "Standard"],
    "AAT": ["N", "Standard"], "AAC": ["N", "Standard"],
    "AAA": ["K", "Standard"], "AAG": ["K", "Standard"],
    "GAT": ["D", "Standard"], "GAC": ["D", "Standard"],
    "GAA": ["E", "Standard"], "GAG": ["E", "Standard"],
    
    "TGT": ["C", "Standard"], "TGC": ["C", "Standard"],
    "TGG": ["W", "Standard"],
    "CGT": ["R", "Standard"], "CGC": ["R", "Standard"], "CGA": ["R", "Standard"], "CGG": ["R", "Standard"],
    "AGT": ["S", "Standard"], "AGC": ["S", "Standard"],
    "AGA": ["S", "Standard"], "AGG": ["S", "Standard"],
    "GGT": ["G", "Standard"], "GGC": ["G", "Standard"], "GGA": ["G", "Standard"], "GGG": ["G", "Standard"]
}
# 6. The Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)
table6 = {
    "TTT": ["F", "Standard"], "TTC": ["F", "Standard"],
    "TTA": ["L", "Standard"], "TTG": ["L", "Standard"],
    "CTT": ["L", "Standard"], "CTC": ["L", "Standard"], "CTA": ["L", "Standard"], "CTG": ["L", "Standard"],
    "ATT": ["I", "Standard"], "ATC": ["I", "Standard"], "ATA": ["I", "Standard"],
    "ATG": ["M", "start"],
    "GTT": ["V", "Standard"], "GTC": ["V", "Standard"], "GTA": ["V", "Standard"], "GTG": ["V", "Standard"],
    
    "TCT": ["S", "Standard"], "TCC": ["S", "Standard"], "TCA": ["S", "Standard"], "TCG": ["S", "Standard"],
    "CCT": ["P", "Standard"], "CCC": ["P", "Standard"], "CCA": ["P", "Standard"], "CCG": ["P", "Standard"],
    "ACT": ["T", "Standard"], "ACC": ["T", "Standard"], "ACA": ["T", "Standard"], "ACG": ["T", "Standard"],
    "GCT": ["A", "Standard"], "GCC": ["A", "Standard"], "GCA": ["A", "Standard"], "GCG": ["A", "Standard"],
    
    "TAT": ["Y", "Standard"], "TAC": ["Y", "Standard"],
    "TAA": ["Q", "Standard"], "TAG": ["Q", "Standard"], "TGA": ["*", "stop"],
    "CAT": ["H", "Standard"], "CAC": ["H", "Standard"],
    "CAA": ["Q", "Standard"], "CAG": ["Q", "Standard"],
    "AAT": ["N", "Standard"], "AAC": ["N", "Standard"],
    "AAA": ["K", "Standard"], "AAG": ["K", "Standard"],
    "GAT": ["D", "Standard"], "GAC": ["D", "Standard"],
    "GAA": ["E", "Standard"], "GAG": ["E", "Standard"],
    
    "TGT": ["C", "Standard"], "TGC": ["C", "Standard"],
    "TGG": ["W", "Standard"],
    "CGT": ["R", "Standard"], "CGC": ["R", "Standard"], "CGA": ["R", "Standard"], "CGG": ["R", "Standard"],
    "AGT": ["S", "Standard"], "AGC": ["S", "Standard"],
    "AGA": ["R", "Standard"], "AGG": ["R", "Standard"],
    "GGT": ["G", "Standard"], "GGC": ["G", "Standard"], "GGA": ["G", "Standard"], "GGG": ["G", "Standard"]
}

# 9. The Echinoderm and Flatworm Mitochondrial Code (transl_table=9)
table9 = {
    "TTT": ["F", "Standard"], "TTC": ["F", "Standard"],
    "TTA": ["L", "Standard"], "TTG": ["L", "start"],
    "CTT": ["L", "Standard"], "CTC": ["L", "Standard"], "CTA": ["L", "Standard"], "CTG": ["L", "Standard"],
    "ATT": ["I", "Standard"], "ATC": ["I", "Standard"], "ATA": ["I", "Standard"],
    "ATG": ["M", "start"],
    "GTT": ["V", "Standard"], "GTC": ["V", "Standard"], "GTA": ["V", "Standard"], "GTG": ["V", "Standard"],
    
    "TCT": ["S", "Standard"], "TCC": ["S", "Standard"], "TCA": ["S", "Standard"], "TCG": ["S", "Standard"],
    "CCT": ["P", "Standard"], "CCC": ["P", "Standard"], "CCA": ["P", "Standard"], "CCG": ["P", "Standard"],
    "ACT": ["T", "Standard"], "ACC": ["T", "Standard"], "ACA": ["T", "Standard"], "ACG": ["T", "Standard"],
    "GCT": ["A", "Standard"], "GCC": ["A", "Standard"], "GCA": ["A", "Standard"], "GCG": ["A", "Standard"],
    
    "TAT": ["Y", "Standard"], "TAC": ["Y", "Standard"],
    "TAA": ["*", "stop"], "TAG": ["*", "stop"], "TGA": ["W", "Standard"],
    "CAT": ["H", "Standard"], "CAC": ["H", "Standard"],
    "CAA": ["N", "Standard"], "CAG": ["N", "Standard"],
    "AAT": ["N", "Standard"], "AAC": ["N", "Standard"],
    "AAA": ["N", "Standard"], "AAG": ["K", "Standard"],
    "GAT": ["D", "Standard"], "GAC": ["D", "Standard"],
    "GAA": ["E", "Standard"], "GAG": ["E", "Standard"],
    
    "TGT": ["C", "Standard"], "TGC": ["C", "Standard"],
    "TGG": ["W", "Standard"],
    "CGT": ["R", "Standard"], "CGC": ["R", "Standard"], "CGA": ["R", "Standard"], "CGG": ["R", "Standard"],
    "AGT": ["S", "Standard"], "AGC": ["S", "Standard"],
    "AGA": ["S", "Standard"], "AGG": ["S", "Standard"],
    "GGT": ["G", "Standard"], "GGC": ["G", "Standard"], "GGA": ["G", "Standard"], "GGG": ["G", "Standard"]
}
# 10. The Euplotid Nuclear Code (transl_table=10)
table10 = {
    "TTT": ["F", "Standard"], "TTC": ["F", "Standard"],
    "TTA": ["L", "Standard"], "TTG": ["L", "Standard"],
    "CTT": ["L", "Standard"], "CTC": ["L", "Standard"], "CTA": ["L", "Standard"], "CTG": ["L", "Standard"],
    "ATT": ["I", "Standard"], "ATC": ["I", "Standard"], "ATA": ["I", "Standard"],
    "ATG": ["M", "start"],
    "GTT": ["V", "Standard"], "GTC": ["V", "Standard"], "GTA": ["V", "Standard"], "GTG": ["V", "Standard"],
    
    "TCT": ["S", "Standard"], "TCC": ["S", "Standard"], "TCA": ["S", "Standard"], "TCG": ["S", "Standard"],
    "CCT": ["P", "Standard"], "CCC": ["P", "Standard"], "CCA": ["P", "Standard"], "CCG": ["P", "Standard"],
    "ACT": ["T", "Standard"], "ACC": ["T", "Standard"], "ACA": ["T", "Standard"], "ACG": ["T", "Standard"],
    "GCT": ["A", "Standard"], "GCC": ["A", "Standard"], "GCA": ["A", "Standard"], "GCG": ["A", "Standard"],
    
    "TAT": ["Y", "Standard"], "TAC": ["Y", "Standard"],
    "TAA": ["*", "stop"], "TAG": ["*", "stop"], "TGA": ["C", "Standard"],
    "CAT": ["H", "Standard"], "CAC": ["H", "Standard"],
    "CAA": ["Q", "Standard"], "CAG": ["Q", "Standard"],
    "AAT": ["N", "Standard"], "AAC": ["N", "Standard"],
    "AAA": ["K", "Standard"], "AAG": ["K", "Standard"],
    "GAT": ["D", "Standard"], "GAC": ["D", "Standard"],
    "GAA": ["E", "Standard"], "GAG": ["E", "Standard"],
    
    "TGT": ["C", "Standard"], "TGC": ["C", "Standard"],
    "TGG": ["W", "Standard"],
    "CGT": ["R", "Standard"], "CGC": ["R", "Standard"], "CGA": ["R", "Standard"], "CGG": ["R", "Standard"],
    "AGT": ["S", "Standard"], "AGC": ["S", "Standard"],
    "AGA": ["R", "Standard"], "AGG": ["R", "Standard"],
    "GGT": ["G", "Standard"], "GGC": ["G", "Standard"], "GGA": ["G", "Standard"], "GGG": ["G", "Standard"]
}

# 11. The Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
table11 = {
    "TTT": ["F", "Standard"], "TCT": ["S", "Standard"], "TAT": ["Y", "Standard"], "TGT": ["C", "Standard"],
    "TTC": ["F", "Standard"], "TCC": ["S", "Standard"], "TAC": ["Y", "Standard"], "TGC": ["C", "Standard"],
    "TTA": ["L", "Standard"], "TCA": ["S", "Standard"], "TAA": ["*", "Stop"], "TGA": ["*", "Stop"],
    "TTG": ["L", "Start"], "TCG": ["S", "Standard"], "TAG": ["*", "Stop"], "TGG": ["W", "Standard"],
    
    "CTT": ["L", "Standard"], "CCT": ["P", "Standard"], "CAT": ["H", "Standard"], "CGT": ["R", "Standard"],
    "CTC": ["L", "Standard"], "CCC": ["P", "Standard"], "CAC": ["H", "Standard"], "CGC": ["R", "Standard"],
    "CTA": ["L", "Standard"], "CCA": ["P", "Standard"], "CAA": ["Q", "Standard"], "CGA": ["R", "Standard"],
    "CTG": ["L", "Start"], "CCG": ["P", "Standard"], "CAG": ["Q", "Standard"], "CGG": ["R", "Standard"],
    
    "ATT": ["I", "Start"], "ACT": ["T", "Standard"], "AAT": ["N", "Standard"], "AGT": ["S", "Standard"],
    "ATC": ["I", "Start"], "ACC": ["T", "Standard"], "AAC": ["N", "Standard"], "AGC": ["S", "Standard"],
    "ATA": ["I", "Start"], "ACA": ["T", "Standard"], "AAA": ["K", "Standard"], "AGA": ["R", "Standard"],
    "ATG": ["M", "Start"], "ACG": ["T", "Standard"], "AAG": ["K", "Standard"], "AGG": ["R", "Standard"],
    
    "GTT": ["V", "Standard"], "GCT": ["A", "Standard"], "GAT": ["D", "Standard"], "GGT": ["G", "Standard"],
    "GTC": ["V", "Standard"], "GCC": ["A", "Standard"], "GAC": ["D", "Standard"], "GGC": ["G", "Standard"],
    "GTA": ["V", "Standard"], "GCA": ["A", "Standard"], "GAA": ["E", "Standard"], "GGA": ["G", "Standard"],
    "GTG": ["V", "Start"], "GCG": ["A", "Standard"], "GAG": ["E", "Standard"], "GGG": ["G", "Standard"]
}
def translateDNA(DataKmer):
    ProtSeq = []
    for i in DataKmer:
        AminoAcid = table11[i]
        if AminoAcid[1] == "Stop":
            ProtSeq = ''.join(ProtSeq)
            return ProtSeq
        else:
            ProtSeq.append(AminoAcid[0])
    ProtSeq = ''.join(ProtSeq)
    return ProtSeq