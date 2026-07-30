#!/usr/bin/env bash

# Script for genomic variant calling pipeline
set -euo pipefail  # Exit on error, undefined variables, and pipe failures

show_help() {
    cat << EOF
Usage: $0 -d <project_dir> -r <ref_fasta> [OPTIONS] <-1 READ1 -2 READ2 | -s SINGLE>

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
  $0 -d ./my_project -r reference.fasta -g annotations.gff -1 reads_R1.fq -2 reads_R2.fq -t 8

  # Single-end reads without annotation
  $0 -d ./output -r ref.fa -s single_reads.fq -p 1 -q 20

EOF
}

# --- Defaults ---
THREADS=4
PLOIDY=2
QUALITY_THRESH=30
PROJECT_DIR=""
REF_FASTA=""
GFF_PATH=""
READ1=""
READ2=""
SINGLE=""

# --- Parse arguments ---
while getopts "d:r:g:t:p:q:1:2:s:h" opt; do
    case $opt in
        d) PROJECT_DIR="$OPTARG" ;;
        r) REF_FASTA="$OPTARG" ;;
        g) GFF_PATH="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        p) PLOIDY="$OPTARG" ;;
        q) QUALITY_THRESH="$OPTARG" ;;
        1) READ1="$OPTARG" ;;
        2) READ2="$OPTARG" ;;
        s) SINGLE="$OPTARG" ;;
        h) show_help; exit 0 ;;
        \?) echo "Error: Invalid option: -$OPTARG" >&2; show_help; exit 1 ;;
        :) echo "Error: Option -$OPTARG requires an argument" >&2; show_help; exit 1 ;;
    esac
done

# --- Validation ---
# Check required parameters
if [[ -z "$PROJECT_DIR" ]]; then
    echo "Error: Project directory (-d) is required" >&2
    show_help
    exit 1
fi

if [[ -z "$REF_FASTA" ]]; then
    echo "Error: Reference FASTA (-r) is required" >&2
    show_help
    exit 1
fi

# Check read input logic
if [[ -n "$SINGLE" && (-n "$READ1" || -n "$READ2") ]]; then
    echo "Error: Cannot specify both single-end (-s) and paired-end (-1/-2) reads" >&2
    exit 1
fi

if [[ -z "$SINGLE" && -z "$READ1" && -z "$READ2" ]]; then
    echo "Error: Must specify either single-end reads (-s) or paired-end reads (-1 and -2)" >&2
    show_help
    exit 1
fi

if [[ -n "$READ1" && -z "$READ2" ]]; then
    echo "Error: -1 (forward reads) requires -2 (reverse reads)" >&2
    exit 1
fi

if [[ -z "$READ1" && -n "$READ2" ]]; then
    echo "Error: -2 (reverse reads) requires -1 (forward reads)" >&2
    exit 1
fi

# Check if files exist
if [[ ! -f "$REF_FASTA" ]]; then
    echo "Error: Reference FASTA file not found: $REF_FASTA" >&2
    exit 1
fi

if [[ -n "$READ1" && ! -f "$READ1" ]]; then
    echo "Error: Forward reads file not found: $READ1" >&2
    exit 1
fi

if [[ -n "$READ2" && ! -f "$READ2" ]]; then
    echo "Error: Reverse reads file not found: $READ2" >&2
    exit 1
fi

if [[ -n "$SINGLE" && ! -f "$SINGLE" ]]; then
    echo "Error: Single-end reads file not found: $SINGLE" >&2
    exit 1
fi

if [[ -n "$GFF_PATH" && ! -f "$GFF_PATH" ]]; then
    echo "Error: GFF annotation file not found: $GFF_PATH" >&2
    exit 1
fi

# Validate numeric parameters
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
    echo "Error: Threads must be a positive integer" >&2
    exit 1
fi

if ! [[ "$PLOIDY" =~ ^[0-9]+$ ]] || [[ "$PLOIDY" -lt 1 ]]; then
    echo "Error: Ploidy must be a positive integer" >&2
    exit 1
fi

if ! [[ "$QUALITY_THRESH" =~ ^[0-9]+$ ]]; then
    echo "Error: Quality threshold must be a non-negative integer" >&2
    exit 1
fi

# --- Setup directory structure ---
REF_PATH=$(realpath "$REF_FASTA")

LOG_DIR="$PROJECT_DIR/logs"
ALIGN_DIR="$PROJECT_DIR/alignments"
SNP_DIR="$PROJECT_DIR/snps"
PROT_DIR="$PROJECT_DIR/proteins"

# Create all required directories
mkdir -p "$PROJECT_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$ALIGN_DIR/Raw"
mkdir -p "$ALIGN_DIR/Filtered"
mkdir -p "$ALIGN_DIR/Indelqual"
mkdir -p "$SNP_DIR/BCFTools"
mkdir -p "$SNP_DIR/Snver"
mkdir -p "$SNP_DIR/Lofreq"

# Create protein directory only if annotation will be performed
if [[ -n "$GFF_PATH" ]]; then
    mkdir -p "$PROT_DIR"
    GFF_PATH=$(realpath "$GFF_PATH")
fi

# --- Display configuration ---
echo "=== Pipeline Configuration ==="
echo "Project directory: $PROJECT_DIR"
echo "Reference FASTA:   $REF_PATH"

if [[ -n "$GFF_PATH" ]]; then
    echo "GFF annotation:    $GFF_PATH"
else
    echo "GFF annotation:    Not provided (skipping annotation)"
fi

if [[ -n "$SINGLE" ]]; then
    echo "Read mode:         Single-end"
    echo "Reads file:        $SINGLE"
else
    echo "Read mode:         Paired-end"
    echo "Forward reads:     $READ1"
    echo "Reverse reads:     $READ2"
fi

echo "Threads:           $THREADS"
echo "Ploidy:            $PLOIDY"
echo "Quality threshold: $QUALITY_THRESH"
echo "=============================="
echo

echo "[INFO] Starting pipeline in: $PROJECT_DIR"
echo "[INFO] Reference: $REF_PATH"
echo "[INFO] Threads: $THREADS | Ploidy: $PLOIDY | Quality ≥ $QUALITY_THRESH"
echo "------------------------------------------------------"

# --- Determine sample name and reads ---
if [[ -n "$SINGLE" ]]; then
    # Single-end mode
    BASENAME=$(basename "$SINGLE")
    SAMPLE=$(echo "$BASENAME" | sed -E 's/\.(fastq|fq)(\.gz)?$//')
    READS=("$SINGLE")
    echo "[INFO] Sample: $SAMPLE → single-end"
else
    # Paired-end mode
    BASENAME=$(basename "$READ1")
    SAMPLE=$(echo "$BASENAME" | sed -E 's/_R1\.(fastq|fq)(\.gz)?$//')
    READS=("$READ1" "$READ2")
    echo "[INFO] Sample: $SAMPLE → paired-end"
fi

# --- Define output paths ---
RAW_SAM="$ALIGN_DIR/Raw/${SAMPLE}.sam"
RAW_BAM="$ALIGN_DIR/Raw/${SAMPLE}.bam"
FILTERED_BAM="$ALIGN_DIR/Filtered/${SAMPLE}.q${QUALITY_THRESH}.bam"
INDELQUAL_BAM="$ALIGN_DIR/Indelqual/${SAMPLE}.indelqual.bam"

# --- Alignment ---
echo "[STEP] Aligning $SAMPLE..."
minimap2 -a -t "$THREADS" "$REF_PATH" "${READS[@]}" > "$RAW_SAM" 2> "$LOG_DIR/${SAMPLE}_minimap2.log"
samtools sort -@ "$THREADS" -o "$RAW_BAM" "$RAW_SAM"
rm "$RAW_SAM"
samtools index "$RAW_BAM"

# --- Quality filter ---
echo "[STEP] Quality Filtering, Quality ≥ $QUALITY_THRESH"
samtools view -q "$QUALITY_THRESH" -b "$RAW_BAM" > "$FILTERED_BAM"
samtools index "$FILTERED_BAM"

# --- LoFreq ---
echo "[STEP] Running LoFreq indel quality..."
lofreq indelqual --dindel --ref "$REF_PATH" --out "$INDELQUAL_BAM" "$RAW_BAM" 2> "$LOG_DIR/${SAMPLE}_indelqual.log"
samtools index "$INDELQUAL_BAM"
# --- LoFreq ---
echo "[STEP] SNP Calling (parallel execution)..."
lofreq call-parallel --pp-threads "$THREADS"\
    -f "$REF_PATH" --call-indels \
    -o "$SNP_DIR/Lofreq/${SAMPLE}.vcf"\
    "$INDELQUAL_BAM" 2> "$LOG_DIR/${SAMPLE}_lofreq.log" &

# --- bcftools ---
bcftools mpileup -Ou -f "$REF_PATH" "$FILTERED_BAM" \
| bcftools call -mv --ploidy "$PLOIDY" -Ov -o "$SNP_DIR/BCFTools/${SAMPLE}.vcf" 2> "$LOG_DIR/${SAMPLE}_bcftools.log" &

# --- snver ---
snver -i "$FILTERED_BAM" -r "$REF_PATH" -o "$SNP_DIR/Snver/${SAMPLE}" -n "$PLOIDY" 2> "$LOG_DIR/${SAMPLE}_snver.log" &
wait
mv "$SNP_DIR/Snver/${SAMPLE}".filter.vcf "$SNP_DIR/Snver/${SAMPLE}".vcf
rm "$SNP_DIR/Snver/${SAMPLE}".failed.log "$SNP_DIR/Snver/${SAMPLE}".indel.filter.vcf "$SNP_DIR/Snver/${SAMPLE}".indel.raw.vcf "$SNP_DIR/Snver/${SAMPLE}".raw.vcf

echo "[DONE] SNP calling completed for $SAMPLE"
echo "------------------------------------------------------"

# --- Normalize VCF files ---
echo "[STEP] Normalizing VCF files..."
TEMP_DIR="$SNP_DIR/Temp"
mkdir -p "$TEMP_DIR"

# Normalize BCFTools VCF
echo "[INFO] Normalizing BCFTools VCF..."
cp "$SNP_DIR/BCFTools/${SAMPLE}.vcf" "$TEMP_DIR/${SAMPLE}.vcf"
bcftools sort "$TEMP_DIR/${SAMPLE}.vcf" -Ov -o "$TEMP_DIR/${SAMPLE}_sort.vcf"
bcftools norm -f "$REF_PATH" -m -both -Ov "$TEMP_DIR/${SAMPLE}_sort.vcf" -o "$TEMP_DIR/${SAMPLE}_norm.vcf"
mv "$TEMP_DIR/${SAMPLE}_norm.vcf" "$SNP_DIR/BCFTools/${SAMPLE}.vcf"
rm "$TEMP_DIR"/*.vcf

# Normalize Strelka VCF
echo "[INFO] Normalizing Snver VCF..."
cp "$SNP_DIR/Snver/${SAMPLE}.vcf" "$TEMP_DIR/${SAMPLE}.vcf"
bgzip "$TEMP_DIR/${SAMPLE}.vcf"
tabix -p vcf "$TEMP_DIR/${SAMPLE}.vcf.gz"
bcftools sort "$TEMP_DIR/${SAMPLE}.vcf.gz" -Ov -o "$TEMP_DIR/${SAMPLE}_sort.vcf"
bcftools norm -f "$REF_PATH" -m -both -Ov "$TEMP_DIR/${SAMPLE}_sort.vcf" -o "$TEMP_DIR/${SAMPLE}_norm.vcf"
mv "$TEMP_DIR/${SAMPLE}_norm.vcf" "$SNP_DIR/Snver/${SAMPLE}.vcf"
rm -f "$TEMP_DIR"/*.vcf "$TEMP_DIR"/*.gz "$TEMP_DIR"/*.tbi

# Normalize Lofreq VCF (requires bgzip and tabix)
echo "[INFO] Normalizing Lofreq VCF..."
cp "$SNP_DIR/Lofreq/${SAMPLE}.vcf" "$TEMP_DIR/${SAMPLE}.vcf"
bgzip "$TEMP_DIR/${SAMPLE}.vcf"
tabix -p vcf "$TEMP_DIR/${SAMPLE}.vcf.gz"
bcftools sort "$TEMP_DIR/${SAMPLE}.vcf.gz" -Ov -o "$TEMP_DIR/${SAMPLE}_sort.vcf"
bcftools norm -f "$REF_PATH" -m -both -Ov "$TEMP_DIR/${SAMPLE}_sort.vcf" -o "$TEMP_DIR/${SAMPLE}_norm.vcf"
mv "$TEMP_DIR/${SAMPLE}_norm.vcf" "$SNP_DIR/Lofreq/${SAMPLE}.vcf"
rm -f "$TEMP_DIR"/*.vcf "$TEMP_DIR"/*.gz "$TEMP_DIR"/*.tbi

# Clean up temp directory
rm -rf "$TEMP_DIR"

echo "[DONE] VCF normalization completed"
echo "------------------------------------------------------"

# --- Annotation (optional) ---
if [[ -n "$GFF_PATH" ]]; then
    echo "[STEP] Running variant annotation..."
    
    # Locate the annotation script
    SCRIPT_DIR="$(dirname "$(realpath "$0")")/Scripts"
    ANNOT_SCRIPT="$SCRIPT_DIR/Annot.py"
    
    if [[ ! -f "$ANNOT_SCRIPT" ]]; then
        echo "[WARNING] Annotation script not found at: $ANNOT_SCRIPT" >&2
        echo "[WARNING] Skipping annotation step" >&2
    else
        python3 "$ANNOT_SCRIPT" \
            -f "$REF_PATH" \
            -g "$GFF_PATH" \
            -v1 "$SNP_DIR/BCFTools/${SAMPLE}.vcf" \
            -v2 "$SNP_DIR/Snver/${SAMPLE}.vcf" \
            -v3 "$SNP_DIR/Lofreq/${SAMPLE}.vcf" \
            -o1 "$PROT_DIR/${SAMPLE}_protein.fasta" \
            -o2 "$PROT_DIR/${SAMPLE}_annotation.tsv" \
            -o3 "$SNP_DIR/${SAMPLE}_consensus.tsv" \
            2> "$LOG_DIR/${SAMPLE}_annot.log"
        
        echo "[DONE] Annotation completed"
    fi
else
    echo "[INFO] Skipping annotation step (no GFF file provided) just merging"
    # Locate the annotation script
    SCRIPT_DIR="$(dirname "$(realpath "$0")")/Scripts"
    ANNOT_SCRIPT="$SCRIPT_DIR/Annot.py"
    
    if [[ ! -f "$ANNOT_SCRIPT" ]]; then
        echo "[WARNING] Annotation script not found at: $ANNOT_SCRIPT" >&2
        echo "[WARNING] Skipping annotation step" >&2
    else
        python3 "$ANNOT_SCRIPT" \
            -f "$REF_PATH" \
            -v1 "$SNP_DIR/BCFTools/${SAMPLE}.vcf" \
            -v2 "$SNP_DIR/Snver/${SAMPLE}.vcf" \
            -v3 "$SNP_DIR/Lofreq/${SAMPLE}.vcf" \
            -o1 "$PROT_DIR/${SAMPLE}_protein.fasta" \
            -o2 "$PROT_DIR/${SAMPLE}_annotation.tsv" \
            -o3 "$SNP_DIR/${SAMPLE}_consensus.tsv" \
            2> "$LOG_DIR/${SAMPLE}_annot.log"
        
        echo "[DONE] Merging completed"
    fi
fi

echo "------------------------------------------------------"
echo "[INFO] Pipeline finished successfully for sample: $SAMPLE"
echo "[INFO] Results are located in: $PROJECT_DIR"
echo "  - Alignments: $ALIGN_DIR"
echo "  - SNPs:       $SNP_DIR"
if [[ -n "$GFF_PATH" ]]; then
    echo "  - Proteins:   $PROT_DIR"
fi
echo "  - Logs:       $LOG_DIR"
