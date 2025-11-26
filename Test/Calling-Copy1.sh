#!/usr/bin/env bash
# ------------------------------------------------------
# Full Variant Calling + Postprocessing Pipeline
# - Uses existing project layout (see below)
# - Optionally accepts a directory with precomputed VCFs (-c)
#
# Expected project layout (exact):
# PROJECT_DIR/
# ├── Input/
# │   ├── Sequencing/
# │   ├── RefGenome/
# │   └── Gff3/
# ├── Output/
# │   ├── Alignments/
# │   │   ├── Raw/
# │   │   ├── Filtered/
# │   │   └── Indelqual/
# │   ├── SNPCalls/
# │   │   ├── BCFTools/
# │   │   ├── Freebayes/
# │   │   └── Lofreq/
# │   └── ProteinSequences/
# │       └── Predictions/
# ├── Logs/
# └── Scripts/
#     └── Annot.py
# ------------------------------------------------------

set -euo pipefail

show_help() {
    cat <<'EOF'
Usage: run_pipeline_full.sh -d <PROJECT_DIR> -r <ref_fasta> -g <gff3> [options]

Required:
  -d <dir>    Project directory (contains Input/, Output/, Logs/, Scripts/)
  -r <file>   Reference FASTA filename (located in Input/RefGenome/)
  -g <file>   GFF3 filename (located in Input/Gff3/)

Options:
  -t <num>    Number of threads (default: 4)
  -p <num>    Ploidy (default: 2)
  -q <num>    MAPQ filter threshold (default: 30)
  -c <dir>    Use existing VCF files from this directory (skip calling)
  -h          Show this help
EOF
}

# defaults
THREADS=4
PLOIDY=2
QUALITY_THRESH=30
VCF_DIR=""

while getopts "d:r:t:p:q:c:g:h" opt; do
    case "$opt" in
        d) PROJECT_DIR="$OPTARG" ;;
        r) REF_FASTA="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        p) PLOIDY="$OPTARG" ;;
        q) QUALITY_THRESH="$OPTARG" ;;
        c) VCF_DIR="$OPTARG" ;;
        g) REF_GFF="$OPTARG" ;;
        h) show_help; exit 0 ;;
        *) show_help; exit 1 ;;
    esac
done

# basic validation
if [[ -z "${PROJECT_DIR:-}" || -z "${REF_FASTA:-}" || -z "${REF_GFF:-}" ]]; then
    echo "[ERROR] Missing required args."
    show_help
    exit 1
fi

# canonical paths
PROJECT_DIR="$(realpath "$PROJECT_DIR")"
REF_PATH="$PROJECT_DIR/Input/RefGenome/$REF_FASTA"
GFF_PATH="$PROJECT_DIR/Input/Gff3/$REF_GFF"
SEQ_DIR="$PROJECT_DIR/Input/Sequencing"
LOG_DIR="$PROJECT_DIR/Logs"
ALIGN_DIR="$PROJECT_DIR/Output/Alignments"
SNP_DIR="$PROJECT_DIR/Output/SNPCalls"
PROT_DIR="$PROJECT_DIR/Output/ProteinSequences/Predictions"
SCRIPTS_DIR="$PROJECT_DIR/Scripts"   # Annot.py is expected here

# ensure target dirs exist (warn if not)
mkdir -p "$LOG_DIR" "$PROT_DIR"
for d in "$SEQ_DIR" "$REF_PATH" "$GFF_PATH" "$ALIGN_DIR/Raw" "$ALIGN_DIR/Filtered" "$ALIGN_DIR/Indelqual" \
         "$SNP_DIR/BCFTools" "$SNP_DIR/Freebayes" "$SNP_DIR/Lofreq" "$SCRIPTS_DIR"; do
    if [[ ! -e "$d" ]]; then
        echo "[WARNING] Expected path missing: $d"
    fi
done

if [[ ! -f "$REF_PATH" ]]; then
    echo "[ERROR] Reference FASTA not found: $REF_PATH"
    exit 1
fi
if [[ ! -f "$GFF_PATH" ]]; then
    echo "[ERROR] GFF3 not found: $GFF_PATH"
    exit 1
fi
if [[ ! -f "$SCRIPTS_DIR/Annot.py" ]]; then
    echo "[ERROR] Annot.py not found in $SCRIPTS_DIR"
    exit 1
fi

echo "[INFO] Project: $PROJECT_DIR"
echo "[INFO] Reference: $REF_PATH"
echo "[INFO] GFF3: $GFF_PATH"
echo "[INFO] Threads: $THREADS | Ploidy: $PLOIDY | MAPQ >= $QUALITY_THRESH"
if [[ -n "$VCF_DIR" ]]; then
    VCF_DIR="$(realpath "$VCF_DIR")"
    echo "[INFO] Using precomputed VCFs from: $VCF_DIR (caller steps will be skipped)"
fi
echo "------------------------------------------------------"

# FASTQ discovery (accepts many Illumina naming variants)
mapfile -t FASTQ_R1_FILES < <(find "$SEQ_DIR" -type f \( -iname "*_R1*.fastq*" -o -iname "*_R1*.fq*" \))
mapfile -t FASTQ_SINGLE < <(find "$SEQ_DIR" -type f \( -iname "*.fastq*" -o -iname "*.fq*" \) -not -iname "*_R1*" -not -iname "*_R2*")

ALL_SAMPLES=( "${FASTQ_R1_FILES[@]}" "${FASTQ_SINGLE[@]}" )
if [[ "${#ALL_SAMPLES[@]}" -eq 0 ]]; then
    echo "[ERROR] No FASTQ files found in $SEQ_DIR"
    exit 1
fi

# helper: find VCF for sample (tries multiple patterns & subdirs)
find_vcfs_for_sample() {
    local sample="$1"
    local vcfdir="$2"
    local bcft="" fb="" lf=""

    # candidate locations: direct in vcfdir, or subdirs BCFTools/ Freebayes/ Lofreq
    # patterns we try: sample*.vcf, *sample*.vcf, sample*.vcf.gz
    for base in "$vcfdir" "$vcfdir/BCFTools" "$vcfdir/Freebayes" "$vcfdir/Lofreq"; do
        [[ -d "$base" ]] || continue
        # bcftools-like
        if [[ -z "$bcft" ]]; then
            bcft_candidate=( "$base/${sample}"*bcftools*.vcf* "$base/${sample}"*.bcf* "$base/${sample}"*.vcf* "$base/"*bcftools*"${sample}"*.vcf* )
            for f in "${bcft_candidate[@]}"; do
                [[ -f "$f" ]] && { bcft="$f"; break; }
            done
        fi
        # freebayes-like
        if [[ -z "$fb" ]]; then
            fb_candidate=( "$base/${sample}"*freebayes*.vcf* "$base/${sample}"*.vcf* "$base/"*freebayes*"${sample}"*.vcf* )
            for f in "${fb_candidate[@]}"; do
                [[ -f "$f" ]] && { fb="$f"; break; }
            done
        fi
        # lofreq-like
        if [[ -z "$lf" ]]; then
            lf_candidate=( "$base/${sample}"*lofreq*.vcf* "$base/${sample}"*.vcf* "$base/"*lofreq*"${sample}"*.vcf* )
            for f in "${lf_candidate[@]}"; do
                [[ -f "$f" ]] && { lf="$f"; break; }
            done
        fi
    done

    # fallback: search anywhere under vcfdir for sample name
    if [[ -z "$bcft" || -z "$fb" || -z "$lf" ]]; then
        while IFS= read -r -d '' f; do
            fname="$(basename "$f")"
            # heuristics: if fname contains 'bcftools' -> bcft; 'freebayes' -> fb; 'lofreq' -> lf
            if [[ -z "$bcft" && "$fname" == *bcftools* ]]; then bcft="$f"; fi
            if [[ -z "$fb" && "$fname" == *freebayes* ]]; then fb="$f"; fi
            if [[ -z "$lf" && "$fname" == *lofreq* ]]; then lf="$f"; fi
            # also accept if filename contains sample name and no caller tag, but prefer exact sample match
            if [[ -z "$bcft" && "$fname" == ${sample}* && "$fname" == *.vcf* ]]; then bcft="$f"; fi
            if [[ -z "$fb" && "$fname" == ${sample}* && "$fname" == *.vcf* ]]; then fb="$f"; fi
            if [[ -z "$lf" && "$fname" == ${sample}* && "$fname" == *.vcf* ]]; then lf="$f"; fi
        done < <(find "$vcfdir" -type f -name "*.vcf*" -print0)
    fi

    # print the three paths (possibly empty) separated by | so caller can parse
    printf "%s|%s|%s" "${bcft:-}" "${fb:-}" "${lf:-}"
}

# main sample loop
for FASTQ_R1 in "${ALL_SAMPLES[@]}"; do
    # normalize path
    FASTQ_R1="$(realpath "$FASTQ_R1")"
    BASENAME="$(basename "$FASTQ_R1")"

    # derive sample name: try to strip _R1*, _R2*, .fastq*, .fq*
    SAMPLE="$(echo "$BASENAME" | sed -E 's/(_R1|_R2)(_.*)?(\.fastq|\.fq)(\.gz)?$//I' | sed -E 's/(\.fastq|\.fq)(\.gz)?$//I')"
    # fallback: remove extensions if previous failed
    SAMPLE="${SAMPLE%%.*}"

    FASTQ_R2="${FASTQ_R1/_R1/_R2}"
    if [[ -f "$FASTQ_R2" ]]; then
        echo "[INFO] $SAMPLE : paired-end detected"
        READS=( "$FASTQ_R1" "$FASTQ_R2" )
    else
        echo "[INFO] $SAMPLE : single-end detected"
        READS=( "$FASTQ_R1" )
    fi

    RAW_SAM="$ALIGN_DIR/Raw/${SAMPLE}.sam"
    RAW_BAM="$ALIGN_DIR/Raw/${SAMPLE}.bam"
    FILTERED_BAM="$ALIGN_DIR/Filtered/${SAMPLE}.q${QUALITY_THRESH}.bam"
    INDELQUAL_BAM="$ALIGN_DIR/Indelqual/${SAMPLE}.indelqual.bam"

    # When VCF_DIR provided -> skip calling and set VCF vars; otherwise run callers
    if [[ -n "$VCF_DIR" ]]; then
        echo "[INFO] Using precomputed VCFs for sample $SAMPLE (skipping calling)"
        # try to find per-sample VCFs inside VCF_DIR
        IFS='|' read -r FOUND_BCFT FOUND_FB FOUND_LF <<< "$(find_vcfs_for_sample "$SAMPLE" "$VCF_DIR")"

        # if not found, try relaxed matching by sample name anywhere in VCF_DIR
        if [[ -z "$FOUND_BCFT" || -z "$FOUND_FB" || -z "$FOUND_LF" ]]; then
            echo "[INFO] Some VCFs missing in VCF_DIR for $SAMPLE, attempting relaxed search..."
            # relaxed search: any file containing sample name
            FOUND_BCFT="$(find "$VCF_DIR" -type f -iname "*${SAMPLE}*.vcf*" | grep -i bcftools -m1 || true)"
            FOUND_FB="$(find "$VCF_DIR" -type f -iname "*${SAMPLE}*.vcf*" | grep -i freebayes -m1 || true)"
            FOUND_LF="$(find "$VCF_DIR" -type f -iname "*${SAMPLE}*.vcf*" | grep -i lofreq -m1 || true)"
            # generic fallback: first file matching sample name for each
            if [[ -z "$FOUND_BCFT" ]]; then FOUND_BCFT="$(find "$VCF_DIR" -type f -iname "*${SAMPLE}*.vcf*" -print -quit || true)"; fi
            if [[ -z "$FOUND_FB"  ]]; then FOUND_FB="$(find "$VCF_DIR" -type f -iname "*${SAMPLE}*.vcf*" -print -quit || true)"; fi
            if [[ -z "$FOUND_LF"  ]]; then FOUND_LF="$(find "$VCF_DIR" -type f -iname "*${SAMPLE}*.vcf*" -print -quit || true)"; fi
        fi

        # final check: at least one VCF must exist to continue with annotation; prefer all three
        if [[ -z "$FOUND_BCFT" && -z "$FOUND_FB" && -z "$FOUND_LF" ]]; then
            echo "[ERROR] No VCFs found in $VCF_DIR for sample $SAMPLE. Skipping."
            continue
        fi

        # choose paths: prioritize found; else use pipeline default locations (if present)
        VCF_BCFT="${FOUND_BCFT:-$SNP_DIR/BCFTools/${SAMPLE}.vcf}"
        VCF_FB="${FOUND_FB:-$SNP_DIR/Freebayes/${SAMPLE}.vcf}"
        VCF_LF="${FOUND_LF:-$SNP_DIR/Lofreq/${SAMPLE}.vcf}"

    else
        # Run alignment & variant calling
        echo "[STEP] Aligning $SAMPLE ..."
        minimap2 -a -t "$THREADS" "$REF_PATH" "${READS[@]}" > "$RAW_SAM" 2> "$LOG_DIR/${SAMPLE}_minimap2.log"
        samtools sort -o "$RAW_BAM" "$RAW_SAM"
        rm -f "$RAW_SAM"
        samtools index "$RAW_BAM"

        echo "[STEP] MAPQ filter (>= $QUALITY_THRESH) ..."
        samtools view -q "$QUALITY_THRESH" -b "$RAW_BAM" > "$FILTERED_BAM"
        samtools index "$FILTERED_BAM"

        echo "[STEP] LoFreq indelqual & call ..."
        lofreq indelqual --dindel --ref "$REF_PATH" --out "$INDELQUAL_BAM" "$RAW_BAM" 2> "$LOG_DIR/${SAMPLE}_indelqual.log"
        samtools index "$INDELQUAL_BAM"

        lofreq call-parallel --pp-threads "$THREADS" -f "$REF_PATH" --call-indels \
            -o "$SNP_DIR/Lofreq/${SAMPLE}.vcf" "$INDELQUAL_BAM" 2> "$LOG_DIR/${SAMPLE}_lofreq.log" &

        echo "[STEP] bcftools mpileup & call ..."
        bcftools mpileup -Ou -f "$REF_PATH" "$FILTERED_BAM" \
            | bcftools call -mv --ploidy "$PLOIDY" -Ov -o "$SNP_DIR/BCFTools/${SAMPLE}.vcf" \
            2> "$LOG_DIR/${SAMPLE}_bcftools.log" &

        echo "[STEP] freebayes ..."
        freebayes -f "$REF_PATH" "$RAW_BAM" --ploidy "$PLOIDY" \
            > "$SNP_DIR/Freebayes/${SAMPLE}.vcf" 2> "$LOG_DIR/${SAMPLE}_freebayes.log" &

        wait
        echo "[INFO] Variant calling finished for $SAMPLE"

        VCF_BCFT="$SNP_DIR/BCFTools/${SAMPLE}.vcf"
        VCF_FB="$SNP_DIR/Freebayes/${SAMPLE}.vcf"
        VCF_LF="$SNP_DIR/Lofreq/${SAMPLE}.vcf"
    fi

    # show which VCFs will be used for annotation
    echo "[INFO] VCFs used for $SAMPLE:"
    echo "       BCFTools: $VCF_BCFT"
    echo "       Freebayes: $VCF_FB"
    echo "       LoFreq: $VCF_LF"

    # require at least one VCF to proceed to annotation
    if [[ ! -f "${VCF_BCFT:-}" && ! -f "${VCF_FB:-}" && ! -f "${VCF_LF:-}" ]]; then
        echo "[WARNING] No VCF files available for $SAMPLE. Skipping annotation."
        continue
    fi

    # Ensure Annot.py finds its modules regardless of calling cwd:
    # We'll call Annot.py from SCRIPTS_DIR and Annot.py should handle sys.path for PythonScripts.
    ANNOT_SCRIPT="$SCRIPTS_DIR/Annot.py"
    if [[ ! -f "$ANNOT_SCRIPT" ]]; then
        echo "[ERROR] Annot.py not found at $ANNOT_SCRIPT. Skipping annotation."
        continue
    fi

    echo "[STEP] Running Annot.py for $SAMPLE ..."
    python3 "$ANNOT_SCRIPT" \
        -f "$REF_PATH" \
        -g "$GFF_PATH" \
        -v1 "${VCF_BCFT:-}" \
        -v2 "${VCF_FB:-}" \
        -v3 "${VCF_LF:-}" \
        -o1 "$PROT_DIR/${SAMPLE}_protein.fasta" \
        -o2 "$PROT_DIR/${SAMPLE}_annotation.tsv" \
        -o3 "$SNP_DIR/${SAMPLE}_consensus.tsv" \
        2> "$LOG_DIR/${SAMPLE}_annot.log" || {
            echo "[ERROR] Annot.py failed for $SAMPLE. Check $LOG_DIR/${SAMPLE}_annot.log"
            continue
        }

    echo "[DONE] Sample $SAMPLE fully processed."
    echo "------------------------------------------------------"
done

echo "[ALL DONE] Pipeline finished for all samples."
