#!/usr/bin/env bash
set -euo pipefail

# Diagnose CRAM issues that cause Manta's "0 high-confidence read pairs" failure.
# Each check maps to a specific failure theory — run against the exact CRAM and
# reference that Manta received.
#
# Usage:
#   ./diagnose_cram_manta.sh <sample.cram> [reference.fasta]
#
# The reference is optional but enables the MD5-mismatch check (theory 5).
# samtools must be on PATH.

CRAM="${1:?Usage: $0 <sample.cram> [reference.fasta]}"
REF="${2:-}"
SAMPLE_N=10000
REGION="chr1:10000000-50000000"

RED='\033[0;31m'
YLW='\033[0;33m'
GRN='\033[0;32m'
BLD='\033[1m'
RST='\033[0m'

pass()  { printf "${GRN}PASS${RST}  %s\n" "$1"; }
warn()  { printf "${YLW}WARN${RST}  %s\n" "$1"; }
fail()  { printf "${RED}FAIL${RST}  %s\n" "$1"; }
header(){ printf "\n${BLD}── Check %s: %s${RST}\n" "$1" "$2"; }

REF_ARG=""
if [[ -n "$REF" ]]; then
    REF_ARG="--reference $REF"
fi

if ! command -v samtools &>/dev/null; then
    echo "Error: samtools not found on PATH" >&2
    exit 1
fi

if [[ ! -f "$CRAM" ]]; then
    echo "Error: CRAM file not found: $CRAM" >&2
    exit 1
fi

printf "${BLD}Diagnosing: %s${RST}\n" "$CRAM"
printf "Sampling %d reads from %s\n" "$SAMPLE_N" "$REGION"

# Grab a sample of reads once; reuse for multiple checks.
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT
SAMPLE_FILE="$TMPDIR/sample.tsv"

samtools view $REF_ARG "$CRAM" "$REGION" 2>"$TMPDIR/view_stderr" \
    | head -n "$SAMPLE_N" > "$SAMPLE_FILE" || true

SAMPLED=$(wc -l < "$SAMPLE_FILE" | tr -d ' ')
if [[ "$SAMPLED" -eq 0 ]]; then
    fail "Could not read any alignments from $REGION — CRAM may be unreadable with this reference or region is empty."
    echo "samtools stderr:"
    cat "$TMPDIR/view_stderr"
    exit 1
fi
printf "Sampled %s reads\n" "$SAMPLED"

# ─────────────────────────────────────────────────────────────────────
header 1 "TLEN (template length) zero or missing"
# If TLEN (col 9) is 0 for all reads, Manta cannot determine pair
# orientation or insert size for any pair.
# ─────────────────────────────────────────────────────────────────────
TLEN_ZERO=$(awk '$9 == 0 {n++} END {print n+0}' "$SAMPLE_FILE")
TLEN_NONZERO=$((SAMPLED - TLEN_ZERO))
TLEN_PCT=$(awk "BEGIN {printf \"%.1f\", ($TLEN_ZERO/$SAMPLED)*100}")

if [[ "$TLEN_NONZERO" -eq 0 ]]; then
    fail "ALL reads have TLEN=0 — Manta cannot determine pair orientation."
elif [[ "$TLEN_ZERO" -gt $((SAMPLED / 2)) ]]; then
    warn "$TLEN_ZERO/$SAMPLED reads ($TLEN_PCT%) have TLEN=0"
else
    pass "$TLEN_NONZERO/$SAMPLED reads have non-zero TLEN ($TLEN_PCT% zero)"
fi

# Show TLEN distribution summary
echo "  TLEN distribution (absolute values, top 5):"
awk '{v=$9; if(v<0) v=-v; print v}' "$SAMPLE_FILE" \
    | sort -n | uniq -c | sort -rn | head -5 \
    | while read count val; do printf "    %8d reads  TLEN=%s\n" "$count" "$val"; done

# ─────────────────────────────────────────────────────────────────────
header 2 "Alignment flags — proper pair, mate unmapped, secondary/supplementary"
# Manta filters out reads that are secondary (0x100), supplementary
# (0x800), have unmapped mates (0x8), or lack the paired flag (0x1).
# If too many reads carry these flags, high-confidence pairs drop to 0.
# ─────────────────────────────────────────────────────────────────────
count_flag() {
    # count reads where (flag & mask) != 0
    local mask=$1
    awk -v m="$mask" '{if(and($2,m)) n++} END {print n+0}' "$SAMPLE_FILE"
}

NOT_PAIRED=$(count_flag 0x1 | xargs -I{} bash -c 'echo $(('"$SAMPLED"' - {}))')
MATE_UNMAP=$(count_flag 8)    # 0x8
SECONDARY=$(count_flag 256)   # 0x100
SUPPLEMENTARY=$(count_flag 2048) # 0x800
PROPER=$(count_flag 2)        # 0x2
DUPLICATE=$(count_flag 1024)  # 0x400

printf "  Paired (0x1):          %d / %d\n" "$((SAMPLED - NOT_PAIRED))" "$SAMPLED"
printf "  Proper pair (0x2):     %d / %d\n" "$PROPER" "$SAMPLED"
printf "  Mate unmapped (0x8):   %d / %d\n" "$MATE_UNMAP" "$SAMPLED"
printf "  Secondary (0x100):     %d / %d\n" "$SECONDARY" "$SAMPLED"
printf "  Supplementary (0x800): %d / %d\n" "$SUPPLEMENTARY" "$SAMPLED"
printf "  Duplicate (0x400):     %d / %d\n" "$DUPLICATE" "$SAMPLED"

EXCLUDED=$((MATE_UNMAP + SECONDARY + SUPPLEMENTARY + DUPLICATE))
USABLE=$((SAMPLED - EXCLUDED))

if [[ "$USABLE" -eq 0 ]]; then
    fail "ALL reads are filtered by flag checks — no usable pairs for Manta."
elif [[ "$PROPER" -eq 0 ]]; then
    fail "ZERO reads have the proper-pair flag (0x2) set."
elif [[ "$USABLE" -lt $((SAMPLED / 2)) ]]; then
    warn "Only $USABLE/$SAMPLED reads survive flag filtering"
else
    pass "$USABLE/$SAMPLED reads survive flag filtering"
fi

# ─────────────────────────────────────────────────────────────────────
header 3 "Mate mapping quality (MQ tag)"
# Manta may require the MQ auxiliary tag. If absent, mate quality is
# treated as 0 and every pair fails the quality filter.
# ─────────────────────────────────────────────────────────────────────
MQ_PRESENT=$(grep -c 'MQ:i:' "$SAMPLE_FILE" || true)
MQ_MISSING=$((SAMPLED - MQ_PRESENT))

if [[ "$MQ_PRESENT" -eq 0 ]]; then
    fail "NO reads carry the MQ (mate mapping quality) tag."
    echo "  Manta may interpret missing MQ as 0, failing all pairs."
elif [[ "$MQ_MISSING" -gt $((SAMPLED / 10)) ]]; then
    warn "$MQ_MISSING/$SAMPLED reads lack the MQ tag"
else
    pass "$MQ_PRESENT/$SAMPLED reads have the MQ tag"
fi

# Also check MC (mate CIGAR) tag — some Manta versions use it
MC_PRESENT=$(grep -c 'MC:Z:' "$SAMPLE_FILE" || true)
printf "  MC (mate CIGAR) tag:   %d / %d reads\n" "$MC_PRESENT" "$SAMPLED"

# ─────────────────────────────────────────────────────────────────────
header 4 "Read group consistency"
# Manta groups reads by RG tag. If header @RG and read RG:Z: tags are
# mismatched or missing, pairing can silently fail.
# ─────────────────────────────────────────────────────────────────────
HEADER_RGS=$(samtools view -H $REF_ARG "$CRAM" 2>/dev/null | grep "^@RG" | grep -oP 'ID:\K[^\t]*' || true)
HEADER_RG_COUNT=$(echo "$HEADER_RGS" | grep -c . || true)

READ_RGS=$(grep -oP 'RG:Z:\K[^\t]*' "$SAMPLE_FILE" | sort -u || true)
READ_RG_COUNT=$(echo "$READ_RGS" | grep -c . || true)
READS_WITH_RG=$(grep -c 'RG:Z:' "$SAMPLE_FILE" || true)
READS_WITHOUT_RG=$((SAMPLED - READS_WITH_RG))

printf "  @RG lines in header:   %d\n" "$HEADER_RG_COUNT"
if [[ -n "$HEADER_RGS" ]]; then
    echo "$HEADER_RGS" | while read rg; do printf "    - %s\n" "$rg"; done
fi
printf "  Distinct RG:Z: on reads: %d\n" "$READ_RG_COUNT"
printf "  Reads missing RG:Z: tag: %d / %d\n" "$READS_WITHOUT_RG" "$SAMPLED"

if [[ "$HEADER_RG_COUNT" -eq 0 && "$READS_WITH_RG" -eq 0 ]]; then
    warn "No read groups in header or on reads (Manta will use empty-string RG)"
elif [[ "$READS_WITHOUT_RG" -eq "$SAMPLED" ]]; then
    fail "Header has @RG but NO reads carry RG:Z: — read group mismatch."
elif [[ "$HEADER_RG_COUNT" -eq 0 && "$READS_WITH_RG" -gt 0 ]]; then
    fail "Reads have RG:Z: tags but header has no @RG lines."
else
    # Check that read RGs are a subset of header RGs
    ORPHAN_RGS=$(comm -23 <(echo "$READ_RGS" | sort) <(echo "$HEADER_RGS" | sort) 2>/dev/null || true)
    if [[ -n "$ORPHAN_RGS" ]]; then
        warn "Some read RG tags not in header: $ORPHAN_RGS"
    else
        pass "Read group tags consistent between header and reads"
    fi
fi

# ─────────────────────────────────────────────────────────────────────
header 5 "CRAM reference MD5 mismatch"
# CRAM encodes reads as diffs against a reference. Decoding with the
# wrong reference silently corrupts base sequences and can cause
# downstream tools to reject reads.
# ─────────────────────────────────────────────────────────────────────
if [[ -z "$REF" ]]; then
    warn "No reference supplied — skipping MD5 check. Re-run with: $0 $CRAM <ref.fasta>"
else
    if [[ ! -f "$REF" ]]; then
        warn "Reference file not found: $REF — skipping"
    else
        # Extract M5 tags from CRAM header
        CRAM_MD5S="$TMPDIR/cram_md5.tsv"
        samtools view -H $REF_ARG "$CRAM" 2>/dev/null \
            | grep "^@SQ" \
            | sed -n 's/.*SN:\([^\t]*\).*M5:\([^\t]*\).*/\1\t\2/p' \
            > "$CRAM_MD5S" || true

        CRAM_SQ_COUNT=$(samtools view -H $REF_ARG "$CRAM" 2>/dev/null | grep -c "^@SQ" || true)
        CRAM_M5_COUNT=$(wc -l < "$CRAM_MD5S" | tr -d ' ')

        if [[ "$CRAM_M5_COUNT" -eq 0 ]]; then
            warn "CRAM header has no M5 tags — cannot verify reference match"
        else
            # Get MD5s from reference .dict or compute from .fai
            REF_DICT="${REF%.fasta}.dict"
            REF_DICT2="${REF%.fa}.dict"
            MISMATCH=0
            CHECKED=0

            if [[ -f "$REF_DICT" || -f "$REF_DICT2" ]]; then
                DICT="${REF_DICT}"
                [[ -f "$DICT" ]] || DICT="$REF_DICT2"

                while IFS=$'\t' read -r sq_name sq_md5; do
                    REF_MD5=$(grep "SN:${sq_name}" "$DICT" 2>/dev/null \
                        | sed -n 's/.*M5:\([a-f0-9]*\).*/\1/p' || true)
                    if [[ -n "$REF_MD5" && "$REF_MD5" != "$sq_md5" ]]; then
                        MISMATCH=$((MISMATCH + 1))
                        if [[ "$MISMATCH" -le 3 ]]; then
                            printf "  Mismatch: %s  CRAM=%s  REF=%s\n" "$sq_name" "$sq_md5" "$REF_MD5"
                        fi
                    fi
                    CHECKED=$((CHECKED + 1))
                done < "$CRAM_MD5S"

                if [[ "$MISMATCH" -gt 0 ]]; then
                    fail "$MISMATCH/$CHECKED contigs have MD5 mismatches — WRONG REFERENCE."
                else
                    pass "All $CHECKED contig MD5s match between CRAM and reference dict"
                fi
            else
                warn "No .dict file found for reference — cannot compare MD5s."
                echo "  Generate one with: samtools dict $REF > ${REF%.fasta}.dict"
            fi
        fi

        # Check contig naming (chr1 vs 1)
        CRAM_CONTIGS=$(samtools view -H $REF_ARG "$CRAM" 2>/dev/null \
            | grep "^@SQ" | head -3 | grep -oP 'SN:\K[^\t]*')
        REF_CONTIGS=$(samtools view -H "$REF" 2>/dev/null \
            | grep "^@SQ" | head -3 | grep -oP 'SN:\K[^\t]*' || \
            head -3 "${REF}.fai" 2>/dev/null | cut -f1 || true)

        if [[ -n "$CRAM_CONTIGS" && -n "$REF_CONTIGS" ]]; then
            CRAM_FIRST=$(echo "$CRAM_CONTIGS" | head -1)
            REF_FIRST=$(echo "$REF_CONTIGS" | head -1)
            if [[ "$CRAM_FIRST" != "$REF_FIRST" ]]; then
                fail "Contig naming mismatch: CRAM has '$CRAM_FIRST', reference has '$REF_FIRST'"
            fi
        fi
    fi
fi

# ─────────────────────────────────────────────────────────────────────
header 6 "CRAM file integrity"
# Check for truncation, missing index, and size sanity.
# ─────────────────────────────────────────────────────────────────────
CRAM_SIZE=$(stat -f%z "$CRAM" 2>/dev/null || stat -c%s "$CRAM" 2>/dev/null || echo 0)
CRAM_SIZE_GB=$(awk "BEGIN {printf \"%.2f\", $CRAM_SIZE / 1073741824}")
printf "  CRAM size: %s bytes (%s GB)\n" "$CRAM_SIZE" "$CRAM_SIZE_GB"

# Check for index
CRAI="${CRAM}.crai"
CRAI_ALT="${CRAM%.cram}.crai"
if [[ -f "$CRAI" ]]; then
    CRAI_SIZE=$(stat -f%z "$CRAI" 2>/dev/null || stat -c%s "$CRAI" 2>/dev/null || echo 0)
    printf "  CRAI index: %s (%s bytes)\n" "$CRAI" "$CRAI_SIZE"
    pass "CRAI index found"
elif [[ -f "$CRAI_ALT" ]]; then
    printf "  CRAI index: %s\n" "$CRAI_ALT"
    pass "CRAI index found (alternate path)"
else
    warn "No .crai index found — Manta requires an index"
fi

# Quick EOF check — try reading the last region
samtools idxstats $REF_ARG "$CRAM" > "$TMPDIR/idxstats.tsv" 2>"$TMPDIR/idx_stderr" && \
    pass "samtools idxstats succeeded (CRAM is structurally intact)" || \
    fail "samtools idxstats failed — CRAM may be truncated or corrupt"

if [[ -s "$TMPDIR/idx_stderr" ]]; then
    echo "  samtools stderr:"
    head -5 "$TMPDIR/idx_stderr" | sed 's/^/    /'
fi

# ─────────────────────────────────────────────────────────────────────
header 7 "DRAGEN signatures and CRAM version"
# DRAGEN CRAMs may use features (quality score binning, proprietary
# compression) that older htslib versions can't fully decode.
# ─────────────────────────────────────────────────────────────────────
HEADER_FULL="$TMPDIR/header.txt"
samtools view -H $REF_ARG "$CRAM" > "$HEADER_FULL" 2>/dev/null || true

DRAGEN_PG=$(grep -i "dragen" "$HEADER_FULL" || true)
if [[ -n "$DRAGEN_PG" ]]; then
    warn "DRAGEN signature found in header — check Manta's htslib compatibility"
    echo "$DRAGEN_PG" | head -3 | sed 's/^/    /'
else
    pass "No DRAGEN signatures in header"
fi

# Check aligner
ALIGNER=$(grep "^@PG" "$HEADER_FULL" | head -5)
if [[ -n "$ALIGNER" ]]; then
    echo "  @PG lines:"
    echo "$ALIGNER" | sed 's/^/    /'
fi

# CRAM version (from first bytes of file)
CRAM_MAGIC=$(xxd -l 6 "$CRAM" 2>/dev/null | head -1 || true)
if echo "$CRAM_MAGIC" | grep -q "4352 414d"; then
    CRAM_VER=$(xxd -s 4 -l 2 -p "$CRAM" 2>/dev/null || true)
    if [[ -n "$CRAM_VER" ]]; then
        MAJOR=$((16#${CRAM_VER:0:2}))
        MINOR=$((16#${CRAM_VER:2:2}))
        printf "  CRAM version: %d.%d\n" "$MAJOR" "$MINOR"
        if [[ "$MAJOR" -ge 3 && "$MINOR" -ge 1 ]]; then
            warn "CRAM 3.1 — older Manta versions may not support this"
        else
            pass "CRAM version compatible with standard htslib"
        fi
    fi
else
    warn "Could not read CRAM magic bytes (might be BAM, not CRAM)"
fi

# ─────────────────────────────────────────────────────────────────────
header 8 "Mate chromosome — inter-chromosomal pairs"
# If all mates map to a different chromosome (col 7 != '='), Manta
# rejects them from orientation stats since it needs same-chromosome
# pairs.
# ─────────────────────────────────────────────────────────────────────
MATE_SAME=$(awk '$7 == "=" {n++} END {print n+0}' "$SAMPLE_FILE")
MATE_DIFF=$(awk '$7 != "=" && $7 != "*" {n++} END {print n+0}' "$SAMPLE_FILE")
MATE_UNMAP_COL=$(awk '$7 == "*" {n++} END {print n+0}' "$SAMPLE_FILE")

printf "  Mate same chrom (=):     %d / %d\n" "$MATE_SAME" "$SAMPLED"
printf "  Mate diff chrom:         %d / %d\n" "$MATE_DIFF" "$SAMPLED"
printf "  Mate unmapped (*):       %d / %d\n" "$MATE_UNMAP_COL" "$SAMPLED"

if [[ "$MATE_SAME" -eq 0 ]]; then
    fail "ALL mates are on a different chromosome — no same-chrom pairs for Manta."
elif [[ "$MATE_SAME" -lt $((SAMPLED / 2)) ]]; then
    warn "Only $MATE_SAME/$SAMPLED mates on same chromosome"
else
    pass "$MATE_SAME/$SAMPLED mates on same chromosome"
fi

# ─────────────────────────────────────────────────────────────────────
header 9 "Duplicate read names (FASTQ merge artifact)"
# If FASTQs from multiple runs were concatenated and runs share read
# names (e.g. same instrument/flowcell ID, or names were simplified),
# the aligner mis-pairs reads. Result: valid flags and MAPQ, but
# nonsensical mate info that Manta rejects.
# ─────────────────────────────────────────────────────────────────────
TOTAL_NAMES=$(awk '{print $1}' "$SAMPLE_FILE" | wc -l | tr -d ' ')
UNIQUE_NAMES=$(awk '{print $1}' "$SAMPLE_FILE" | sort -u | wc -l | tr -d ' ')
DUP_NAMES=$((TOTAL_NAMES - UNIQUE_NAMES))

# For paired reads, each name appears exactly twice (R1+R2 at same locus).
# More than 2 occurrences means name collisions from merged runs.
NAMES_GT2=$(awk '{print $1}' "$SAMPLE_FILE" | sort | uniq -c | awk '$1 > 2 {n++} END {print n+0}')

printf "  Total read names:      %d\n" "$TOTAL_NAMES"
printf "  Unique read names:     %d\n" "$UNIQUE_NAMES"
printf "  Names appearing >2x:   %d\n" "$NAMES_GT2"

if [[ "$NAMES_GT2" -gt 0 ]]; then
    fail "$NAMES_GT2 read names appear more than twice — likely FASTQ merge collision."
    echo "  Examples of colliding names:"
    awk '{print $1}' "$SAMPLE_FILE" | sort | uniq -c | sort -rn \
        | awk '$1 > 2' | head -5 \
        | while read count name; do printf "    %s  (%dx)\n" "$name" "$count"; done

    # Show whether colliding reads map to different positions
    echo "  Positions for top colliding name:"
    TOP_COLLISION=$(awk '{print $1}' "$SAMPLE_FILE" | sort | uniq -c | sort -rn | awk '$1>2 {print $2; exit}')
    awk -v n="$TOP_COLLISION" '$1 == n {printf "    flag=%s chr=%s pos=%s mate_chr=%s mate_pos=%s tlen=%s\n", $2, $3, $4, $7, $8, $9}' "$SAMPLE_FILE"
else
    pass "No read name collisions detected in sample"
fi

# Also check read name structure — are names from multiple instruments/runs?
echo "  Read name prefixes (instrument:run:flowcell — expect 1 if single run):"
awk '{split($1,a,":"); if(length(a)>=3) print a[1]":"a[2]":"a[3]; else print $1}' "$SAMPLE_FILE" \
    | sort | uniq -c | sort -rn | head -5 \
    | while read count prefix; do printf "    %8d reads  %s\n" "$count" "$prefix"; done

PREFIX_COUNT=$(awk '{split($1,a,":"); if(length(a)>=3) print a[1]":"a[2]":"a[3]; else print $1}' "$SAMPLE_FILE" \
    | sort -u | wc -l | tr -d ' ')
if [[ "$PREFIX_COUNT" -gt 1 ]]; then
    warn "$PREFIX_COUNT distinct instrument:run:flowcell prefixes — data was merged from multiple runs"
else
    pass "Single sequencing run detected"
fi

# ─────────────────────────────────────────────────────────────────────
header 10 "Insert size distribution shape (merged library detection)"
# If FASTQs from libraries with different insert sizes were merged,
# the distribution is multimodal. Manta expects a unimodal distribution
# to determine pair orientation.
# ─────────────────────────────────────────────────────────────────────
# Compute insert size stats from same-chrom, properly-paired reads
ISIZE_FILE="$TMPDIR/isizes.txt"
awk '$7 == "=" && and($2,0x2) && $9 > 0 {print $9}' "$SAMPLE_FILE" > "$ISIZE_FILE"
ISIZE_COUNT=$(wc -l < "$ISIZE_FILE" | tr -d ' ')

if [[ "$ISIZE_COUNT" -lt 100 ]]; then
    warn "Only $ISIZE_COUNT positive insert sizes from proper pairs — too few to assess distribution"
else
    ISIZE_MEAN=$(awk '{s+=$1} END {printf "%.0f", s/NR}' "$ISIZE_FILE")
    ISIZE_MEDIAN=$(sort -n "$ISIZE_FILE" | awk -v n="$ISIZE_COUNT" 'NR==int(n/2) {print}')
    ISIZE_STDDEV=$(awk -v m="$ISIZE_MEAN" '{d=$1-m; s+=d*d} END {printf "%.0f", sqrt(s/NR)}' "$ISIZE_FILE")
    ISIZE_MIN=$(head -1 <(sort -n "$ISIZE_FILE"))
    ISIZE_MAX=$(tail -1 <(sort -n "$ISIZE_FILE"))
    ISIZE_CV=$(awk "BEGIN {printf \"%.2f\", $ISIZE_STDDEV / $ISIZE_MEAN}")

    printf "  Proper-pair insert sizes (n=%d):\n" "$ISIZE_COUNT"
    printf "    Mean:    %s\n" "$ISIZE_MEAN"
    printf "    Median:  %s\n" "$ISIZE_MEDIAN"
    printf "    Std dev: %s\n" "$ISIZE_STDDEV"
    printf "    Min:     %s\n" "$ISIZE_MIN"
    printf "    Max:     %s\n" "$ISIZE_MAX"
    printf "    CV:      %s\n" "$ISIZE_CV"

    # Check for multimodality: bin into 50bp windows, find peaks
    echo "  Insert size histogram (50bp bins):"
    awk '{bin=int($1/50)*50; bins[bin]++} END {for(b in bins) print b, bins[b]}' "$ISIZE_FILE" \
        | sort -n | awk -v max=0 '{if($2>max) max=$2} END {CONVFMT="%.0f"} {printf "    %4d-%4dbp: %6d %s\n", $1, $1+49, $2, ""}' \
        > "$TMPDIR/hist_raw.txt" || true

    # Re-do with proper max calculation for bar chart
    awk '{bin=int($1/50)*50; bins[bin]++} END {for(b in bins) print b, bins[b]}' "$ISIZE_FILE" \
        | sort -n > "$TMPDIR/hist_data.txt"

    HIST_MAX=$(awk '{if($2>m)m=$2} END{print m}' "$TMPDIR/hist_data.txt")
    awk -v mx="$HIST_MAX" '{
        barlen = int($2 / mx * 40)
        bar = ""
        for(i=0; i<barlen; i++) bar = bar "#"
        printf "    %4d-%4dbp: %6d  %s\n", $1, $1+49, $2, bar
    }' "$TMPDIR/hist_data.txt"

    # Detect multimodality: count bins that are local maxima
    PEAKS=$(awk '
        NR>1 {prev2=prev1; prev1=$2}
        NR>2 {if(prev1 > prev2 && prev1 > $2 && prev1 > 0.1*mx) peaks++}
        {if($2>mx) mx=$2}
        END {print peaks+0}
    ' "$TMPDIR/hist_data.txt")

    if [[ "$PEAKS" -gt 1 ]]; then
        fail "Insert size distribution has $PEAKS peaks — likely merged libraries with different insert sizes."
    elif (( $(awk "BEGIN {print ($ISIZE_CV > 0.40) ? 1 : 0}") )); then
        warn "High coefficient of variation ($ISIZE_CV) — insert sizes are unusually spread"
    else
        pass "Insert size distribution appears unimodal (CV=$ISIZE_CV)"
    fi
fi

# ─────────────────────────────────────────────────────────────────────
header 11 "Pre-delivery deduplication (suspiciously low dup rate)"
# If duplicates were physically removed from FASTQs before alignment,
# the dup rate will be impossibly low for the coverage depth. This
# hints at upstream processing that may have also broken pair info.
# ─────────────────────────────────────────────────────────────────────
DUP_FLAGGED=$(count_flag 1024)
DUP_PCT=$(awk "BEGIN {printf \"%.3f\", ($DUP_FLAGGED / $SAMPLED) * 100}")

printf "  Duplicate-flagged reads: %d / %d (%s%%)\n" "$DUP_FLAGGED" "$SAMPLED" "$DUP_PCT"

# Estimate coverage from the sampled region to gauge expected dup rate.
# At 30-50x WGS, <1% duplication is physically implausible for a single
# PCR library and very unlikely even for PCR-free.
if (( $(awk "BEGIN {print ($DUP_PCT < 1.0) ? 1 : 0}") )); then
    warn "Duplication rate is ${DUP_PCT}% — if this is >=30x WGS, duplicates were likely removed upstream."
    echo "  Expected: 5-15% for PCR-free, 15-30% for PCR libraries at 30-50x."
    echo "  Pre-delivery dedup can break mate consistency or orphan reads."
else
    pass "Duplication rate ($DUP_PCT%) is in plausible range"
fi

# Check for orphaned reads — paired flag set but mate info is missing/inconsistent
PAIRED_COUNT=$(count_flag 1)
ORPHAN_CANDIDATES=$(awk 'and($2,0x1) && ($7 == "*" || $8 == 0) && !and($2,0x4) && !and($2,0x8) {n++} END {print n+0}' "$SAMPLE_FILE")
if [[ "$ORPHAN_CANDIDATES" -gt 0 ]]; then
    ORPHAN_PCT=$(awk "BEGIN {printf \"%.2f\", ($ORPHAN_CANDIDATES / $SAMPLED) * 100}")
    warn "$ORPHAN_CANDIDATES reads ($ORPHAN_PCT%) are flagged as paired but have missing/zero mate position"
    echo "  These may be orphans from inconsistent upstream deduplication."
else
    pass "No orphaned-pair signatures detected"
fi

# ─────────────────────────────────────────────────────────────────────
printf "\n${BLD}── Summary${RST}\n"
printf "If all checks pass, the issue may be Manta-version-specific.\n"
printf "Try converting to BAM first:  samtools view -b -o sample.bam %s\n" "$CRAM"
printf "Then re-run Manta with the BAM to isolate CRAM-specific issues.\n"