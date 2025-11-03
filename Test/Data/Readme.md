# Project Data Directory

This is your working directory.

## Structure
- **Input/** → Store raw input files (e.g., raw reads, reference genome, annotations)
  - Gff3/ → GFF3 annotation files
  - RefGenome/ → Reference genome FASTA files
  - Sequencing/ → Raw sequencing reads
- **Output/** → Contains results from analysis tools
  - Alignments/ → Sequence alignment results
  - ProteinSequences/ → Protein sequence results
  - SNPCalls/ → SNP calling results from various tools
  - Logs/ → Logging data
\n Directory structure:
./Data
├── Input
│   ├── Gff3
│   ├── RefGenome
│   └── Sequencing
├── Output
│   ├── Alignments
│   │   ├── Filtered
│   │   ├── Indelqual
│   │   └── Raw
│   ├── Logs
│   ├── ProteinSequences
│   │   ├── BLAST_Hits
│   │   └── Predictions
│   └── SNPCalls
│       ├── BCFTools
│       ├── Freebayes
│       ├── Lofreq
│       └── Snver
└── Readme.md

19 directories, 1 file
