#!/usr/bin/awk -f

BEGIN {
    FS = "\t"
    OFS = "\t"
    print "POS\tCHROM\tREF\tALT\tSTATUS"
}

# Kommentarzeilen überspringen
/^#/ { next }

{
    chrom  = $1
    pos    = $2
    ref    = $4
    alt    = $5

    # STATUS bestimmen: INDEL wenn REF oder ALT länger als 1 Zeichen
    if (length(ref) > 1 || length(alt) > 1)
        status = "INDEL"
    else
        status = "SNP"

    print pos, chrom, ref, alt, status
}
