#!/bin/bash

# ------------------------------
# Setup Bioinformatics Project Directory
# Usage: ./setup_directory.sh -d /path/to/project
# ------------------------------

# Parse arguments
while getopts "d:" opt; do
    case $opt in
        d) PROJECT_DIR="$OPTARG" ;;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        :) echo "Option -$OPTARG requires a directory path." >&2; exit 1 ;;
    esac
done

# Check if directory was provided
if [ -z "$PROJECT_DIR" ]; then
    echo "Usage: $0 -d /path/to/project"
    exit 1
fi

echo "Setting up the working directory at: $PROJECT_DIR"

# Ensure the project directory exists
mkdir -p "$PROJECT_DIR"
cd "$PROJECT_DIR" || { echo "Error: Cannot access $PROJECT_DIR"; exit 1; }

# Create Data folder and subdirectories
mkdir -p Data/Input/{Gff3,RefGenome,Sequencing}
mkdir -p Data/Output/Alignments/{Filtered,Indelqual,Raw}
mkdir -p Data/Output/ProteinSequences/{BLAST_Hits,Predictions}
mkdir -p Data/Output/SNPCalls/{BCFTools,Freebayes,Lofreq,Snver}
mkdir -p Data/Output/Logs

# Create Readme.md with description
README="$PROJECT_DIR/Data/Readme.md"
{
    echo "# Project Data Directory"
    echo
    echo "This is your working directory."
    echo
    echo "## Structure"
    echo "- **Input/** → Store raw input files (e.g., raw reads, reference genome, annotations)"
    echo "  - Gff3/ → GFF3 annotation files"
    echo "  - RefGenome/ → Reference genome FASTA files"
    echo "  - Sequencing/ → Raw sequencing reads"
    echo "- **Output/** → Contains results from analysis tools"
    echo "  - Alignments/ → Sequence alignment results"
    echo "  - ProteinSequences/ → Protein sequence results"
    echo "  - SNPCalls/ → SNP calling results from various tools"
    echo "  - Logs/ → Logging data"
} > "$README"

echo "Setup complete!"
echo "Readme created at: $README"

# Optional: show folder tree if 'tree' is installed
if command -v tree &> /dev/null; then
    echo "\n Directory structure:" >> $README
    tree "$PROJECT_DIR/Data" >> $README
fi
