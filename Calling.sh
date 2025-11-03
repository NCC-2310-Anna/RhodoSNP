#!/bin/bash
# ------------------------------------------------------
# Variant Calling Pipeline (Paired/Single-End)
# Compatible with existing project structure:
# Input/{Sequencing,RefGenome,Gff3}
# Output/{Alignments,SNPCalls,ProteinSequences}
# ------------------------------------------------------

show_help() {
    echo "Usage: $0 [-d <project_dir>] [-r <ref_fasta>] [-t <threads>] [-p <ploidy>] [-q <quality>] [-h]"
    echo
    echo "Options:"
    echo "  -d <dir>    Base project directory (containing Input/, Output/, Logs/)"
    echo "  -r <ref>    Reference FASTA filename (in Input/RefGenome/)"
    echo "  -t <num>    Number of threads (default: 4) for Lofreq"
    echo "  -p <num>    Ploidy (default: 2)"
    echo "  -q <num>    Mapping quality threshold (default: 30)"
    echo "  -g <gff>    GFF3 filename (in Input/Gff3/)"
    echo "  -h          Show this help message"
}

# --- Defaults ---
THREADS=4
PLOIDY=2
QUALITY_THRESH=30

while getopts "d:r:t:p:q:g:h" opt; do
    case $opt in
        d) PROJECT_DIR="$OPTARG" ;;
        r) REF_FASTA="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        p) PLOIDY="$OPTARG" ;;
        q) QUALITY_THRESH="$OPTARG" ;;
        g) REF_GFF="$OPTARG" ;;
        h) show_help; exit 0 ;;
        *) show_help; exit 1 ;;
    esac
done

# --- Checks ---
if [[ -z "$PROJECT_DIR" || -z "$REF_FASTA" || -z "$REF_GFF" ]]; then
    echo "[ERROR] Missing required arguments."; show_help; exit 1
fi

REF_PATH="$PROJECT_DIR/Input/RefGenome/$REF_FASTA"
SEQ_DIR="$PROJECT_DIR/Input/Sequencing"
LOG_DIR="$PROJECT_DIR/Logs"
ALIGN_DIR="$PROJECT_DIR/Output/Alignments"
SNP_DIR="$PROJECT_DIR/Output/SNPCalls"
GFF_PATH="$PROJECT_DIR/Input/Gff3/$REF_GFF"

if [[ ! -f "$REF_PATH" ]]; then
    echo "[ERROR] Reference FASTA not found: $REF_PATH"; exit 1
fi

if [[ ! -f "$GFF_PATH" ]]; then
    echo "[ERROR] Reference GFF not found: $GFF_PATH"; exit 1
fi

cd "$PROJECT_DIR" || exit 1

# --- Verify directories ---
REQUIRED_DIRS=("$SEQ_DIR" "$LOG_DIR" "$ALIGN_DIR/Raw" "$ALIGN_DIR/Filtered" "$ALIGN_DIR/Indelqual" \
               "$SNP_DIR/BCFTools" "$SNP_DIR/Freebayes" "$SNP_DIR/Lofreq")

for d in "${REQUIRED_DIRS[@]}"; do
    [[ ! -d "$d" ]] && echo "[WARNING] Directory missing: $d"
done

echo "[INFO] Starting pipeline in: $PROJECT_DIR"
echo "[INFO] Reference: $REF_PATH"
echo "[INFO] Gff3: $GFF_PATH"
echo "[INFO] Threads: $THREADS | Ploidy: $PLOIDY | Quality ≥ $QUALITY_THRESH"
echo "------------------------------------------------------"

# --- FASTQ detection ---
FASTQ_R1_FILES=($(find "$SEQ_DIR" -type f \( -name "*_R1.fq*" -o -name "*_R1.fastq*" \)))
echo $FASTQ_R1_FILES
FASTQ_SINGLE=($(find "$SEQ_DIR" -type f \( -regex ".*\.fastq(\.gz)?" \) ! -regex ".*_R[12](_[0-9]+)?\.fastq(\.gz)?" ))

ALL_SAMPLES=("${FASTQ_R1_FILES[@]}" "${FASTQ_SINGLE[@]}")
[[ ${#ALL_SAMPLES[@]} -eq 0 ]] && { echo "[ERROR] No FASTQ files found in $SEQ_DIR"; exit 1; }

# --- Process each sample ---
for FASTQ_R1 in "${ALL_SAMPLES[@]}"; do
    BASENAME=$(basename "$FASTQ_R1")
    SAMPLE=$(echo "$BASENAME" | sed -E 's/_R1|_R2|.fastq(.gz)?//g')
    FASTQ_R2="${FASTQ_R1/_R1/_R2}"

    if [[ -f "$FASTQ_R2" ]]; then
        echo "[INFO] Sample: $SAMPLE → paired-end"
        READS=("$FASTQ_R1" "$FASTQ_R2")
    else
        echo "[INFO] Sample: $SAMPLE → single-end"
        READS=("$FASTQ_R1")
    fi

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
    samtools view -q "$QUALITY_THRESH" -b "$RAW_BAM" > "$FILTERED_BAM"
    samtools index "$FILTERED_BAM"

    # --- LoFreq ---
    lofreq indelqual --dindel --ref "$REF_PATH" --out "$INDELQUAL_BAM" "$RAW_BAM" 2> "$LOG_DIR/${SAMPLE}_indelqual.log"
    samtools index "$INDELQUAL_BAM"

    lofreq call-parallel --pp-threads "$THREADS_PER_CALLER" \
        -f "$REF_PATH" --call-indels \
        -o "$SNP_DIR/Lofreq/${SAMPLE}.vcf" \
        "$INDELQUAL_BAM" 2> "$LOG_DIR/${SAMPLE}_lofreq.log" &

    # --- bcftools ---
    bcftools mpileup -Ou -f "$REF_PATH" "$FILTERED_BAM" \
    | bcftools call -mv --ploidy "$PLOIDY" -Ov -o "$SNP_DIR/BCFTools/${SAMPLE}.vcf" 2> "$LOG_DIR/${SAMPLE}_bcftools.log" &

    # --- FreeBayes ---
    freebayes -f "$REF_PATH" "$RAW_BAM" --ploidy "$PLOIDY" \
        > "$SNP_DIR/Freebayes/${SAMPLE}.vcf" 2> "$LOG_DIR/${SAMPLE}_freebayes.log" &

    wait
    echo "[DONE] $SAMPLE completed."
    echo "------------------------------------------------------"
done

echo "[INFO] Pipeline finished for all samples."
echo "[STEP] Translating SNPs to protein changes..."
python3 $PROJECT_DIR/../Scripts/Annot.py \
  -f "$PROJECT_DIR/Input/RefGenome/$REF_FASTA" \
  -g "$GFF_PATH" \
  -v1 "$SNP_DIR/BCFTools/${SAMPLE}.vcf" \
  -v2 "$SNP_DIR/Freebayes/${SAMPLE}.vcf" \
  -v3 "$SNP_DIR/Lofreq/LF_${SAMPLE}.vcf" \
  -o1 "$PROJECT_DIR/Output/ProteinSequences/Predictions/${SAMPLE}_protein.fasta" \
  -o2 "$PROJECT_DIR/Output/ProteinSequences/Predictions/${SAMPLE}_annotation.tsv" \
  -o3 "$PROJECT_DIR/Output/SNPCalls/${SAMPLE}_consensus.tsv" \
  2> "$LOG_DIR/${SAMPLE}_variant_to_protein.log"


