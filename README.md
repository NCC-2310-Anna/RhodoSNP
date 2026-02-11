# RhodoSNP
An all-in-one variant calling and annotation pipeline for Rhodobacter sphaeroides

## Getting started

```bash
git clone https://github.com/NCC-2310-Anna/RhodoSNP.git
cd RhodoSNP/Setup
mamba env create --name <env name> --file=RhodoSNP_Setup.yml # or conda but that is slower
#make RhodoSNP executable
cd RhodoSNP/Pipeline
chmod +x Calling.sh
```
The Output should look something like that
```
Usage: ./Pipeline/Calling.sh -d <project_dir> -r <ref_fasta> [OPTIONS] <-1 READ1 -2 READ2 | -s SINGLE>

A genomics variant calling pipeline wrapper.

Required Options:
  -d <dir>    Output project directory
  -r <ref>    Reference FASTA path
  
  Read Input (one of the following):
  -1 <file>   Forward reads (paired-end, requires -2)
  -2 <file>   Reverse reads (paired-end, requires -1)
  -s <file>   Single-end reads

Optional Parameters:
  -g <file>   GFF annotation file (enables annotation step)
  -t <num>    Number of threads (default: 4)
  -p <num>    Ploidy (default: 2)
  -q <num>    Mapping quality threshold (default: 30)
  -h          Show this help message

Examples:
  # Paired-end reads with annotation
  ./Pipeline/Calling.sh -d ./my_project -r reference.fasta -g annotations.gff -1 reads_R1.fq -2 reads_R2.fq -t 8

  # Single-end reads without annotation
  ./Pipeline/Calling.sh -d ./output -r ref.fa -s single_reads.fq -p 1 -q 20
```
You can now run the pipeline with the following command
```bash
path/to/RhodoSNP/Pipeline/Calling.sh -h
```

## Use cases

### SNP Calling without further annotation
 ```bash
./Pipeline/Calling.sh -d ./output -r ref.fa -s single_reads.fq -p 1 -q 20
```

### SNP Calling with additional SNP annotation
 ```bash
./Pipeline/Calling.sh -d ./my_project -r reference.fasta -g annotations.gff -1 reads_R1.fq -2 reads_R2.fq -t 8
```
