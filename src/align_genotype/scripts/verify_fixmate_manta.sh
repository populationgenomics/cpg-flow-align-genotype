#!/usr/bin/env bash
set -eu

# Verify what is causing Manta's "0 high-confidence read pairs" failure.
#
# 1. Checks reference compatibility (masked vs unmasked)
# 2. Extracts a region to BAM and attempts Manta on it
# 3. Reports whether the issue is reference, CRAM format, or something else
#
# Usage:
#   ./verify_fixmate_manta.sh <sample.cram> <reference.fasta> [output_bam_path]

CRAM="${1:?Usage: $0 <sample.cram> <reference.fasta> [output_bam_path]}"
REF="${2:?Usage: $0 <sample.cram> <reference.fasta> [output_bam_path]}"
OUTPUT_BAM="${3:-}"
REGION="chr1:10000000-50000000"
SAMPLE_N=10000
CPU=3

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

# ─────────────────────────────────────────────────────────────────────
echo "=== Check 1: Reference compatibility ==="
# ─────────────────────────────────────────────────────────────────────
echo "Reference provided: $REF"

# What reference was the CRAM encoded against?
CRAM_REF=$(samtools view -H --reference "$REF" "$CRAM" 2>/dev/null \
    | grep "^@PG" | grep -oP 'CL:samtools view.*-T \K[^ ]+' || true)
if [ -n "$CRAM_REF" ]; then
    echo "CRAM was encoded with reference: $CRAM_REF"
else
    echo "Could not determine encoding reference from @PG headers"
fi

# Compare first few contig MD5s between CRAM header and provided reference
echo ""
echo "Comparing contig MD5s (CRAM header vs provided reference)..."
CRAM_MD5=$(samtools view -H --reference "$REF" "$CRAM" 2>/dev/null \
    | grep "^@SQ" | grep "SN:chr1\b" \
    | sed -n 's/.*M5:\([a-f0-9]*\).*/\1/p' || true)

if [ -f "${REF%.fasta}.dict" ]; then
    REF_DICT="${REF%.fasta}.dict"
elif [ -f "${REF%.fa}.dict" ]; then
    REF_DICT="${REF%.fa}.dict"
else
    REF_DICT=""
fi

if [ -n "$REF_DICT" ] && [ -n "$CRAM_MD5" ]; then
    REF_MD5=$(grep -P "SN:chr1\t" "$REF_DICT" 2>/dev/null \
        | sed -n 's/.*M5:\([a-f0-9]*\).*/\1/p' || true)
    echo "  chr1 CRAM MD5: $CRAM_MD5"
    echo "  chr1 REF  MD5: $REF_MD5"
    if [ "$CRAM_MD5" = "$REF_MD5" ]; then
        echo "  MATCH — provided reference is compatible with CRAM"
    else
        echo "  MISMATCH — provided reference does NOT match CRAM encoding reference!"
        echo "  This would cause Manta to fail silently."
    fi
else
    echo "  Could not compare (no .dict file or no M5 tags in CRAM header)"
fi

# Check if reference is masked or unmasked
echo ""
REF_BASENAME=$(basename "$REF")
echo "Reference filename: $REF_BASENAME"
if echo "$REF_BASENAME" | grep -qi "masked"; then
    echo "  Reference appears to be MASKED (matches DRAGEN convention)"
else
    echo "  WARNING: Reference does not have 'masked' in filename."
    echo "  CRAM was encoded with Homo_sapiens_assembly38_masked.fasta (from @PG header)."
    echo "  If Manta uses an unmasked reference, that is the likely cause of failure."
fi

# ─────────────────────────────────────────────────────────────────────
echo ""
echo "=== Check 2: Extract region to BAM ==="
# ─────────────────────────────────────────────────────────────────────
samtools view -b --reference "$REF" -o "$TMPDIR/region.bam" "$CRAM" "$REGION"
samtools index "$TMPDIR/region.bam"
TOTAL=$(samtools view -c "$TMPDIR/region.bam")
echo "Extracted $TOTAL reads from $REGION to BAM"

# Verify the BAM is readable without a reference (BAM is self-contained)
BAM_CHECK=$(samtools view -c "$TMPDIR/region.bam" 2>/dev/null || echo "FAIL")
echo "BAM self-contained read check: $BAM_CHECK reads"

# ─────────────────────────────────────────────────────────────────────
echo ""
echo "=== Check 3: Attempt Manta on extracted BAM ==="
# ─────────────────────────────────────────────────────────────────────
MANTA_BIN=""
if [ -f /usr/local/share/manta/bin/configManta.py ]; then
    MANTA_BIN=/usr/local/share/manta/bin/configManta.py
elif command -v configManta.py &>/dev/null; then
    MANTA_BIN=$(command -v configManta.py)
fi

if [ -n "$MANTA_BIN" ]; then
    echo "Found Manta at: $MANTA_BIN"

    echo ""
    echo "--- 3a: Manta on BAM (bypasses CRAM decoding) ---"
    mkdir -p "$TMPDIR/manta_bam"
    python "$MANTA_BIN" \
        --bam "$TMPDIR/region.bam" \
        --referenceFasta "$REF" \
        --runDir "$TMPDIR/manta_bam" 2>&1 || true

    if [ -f "$TMPDIR/manta_bam/runWorkflow.py" ]; then
        python "$TMPDIR/manta_bam/runWorkflow.py" -j "$CPU" 2>&1 || true
        BAM_EXIT=$?
        echo "Manta on BAM exit code: $BAM_EXIT"
        if [ -f "$TMPDIR/manta_bam/workspace/pyflow.data/logs/pyflow_log.txt" ]; then
            echo "--- Manta BAM log (last 30 lines) ---"
            tail -30 "$TMPDIR/manta_bam/workspace/pyflow.data/logs/pyflow_log.txt"
        fi
        # Check for the specific error
        if grep -q "0 high-confidence" "$TMPDIR/manta_bam/workspace/pyflow.data/logs/pyflow_log.txt" 2>/dev/null; then
            echo "RESULT: Manta STILL fails on BAM — issue is NOT CRAM-specific"
        else
            echo "RESULT: Manta succeeded on BAM — issue IS CRAM-specific"
        fi
    else
        echo "configManta failed — check output above"
    fi

    echo ""
    echo "--- 3b: Manta on CRAM directly (for comparison) ---"
    mkdir -p "$TMPDIR/manta_cram"
    python "$MANTA_BIN" \
        --bam "$CRAM" \
        --referenceFasta "$REF" \
        --runDir "$TMPDIR/manta_cram" 2>&1 || true

    if [ -f "$TMPDIR/manta_cram/runWorkflow.py" ]; then
        python "$TMPDIR/manta_cram/runWorkflow.py" -j "$CPU" 2>&1 || true
        CRAM_EXIT=$?
        echo "Manta on CRAM exit code: $CRAM_EXIT"
        if [ -f "$TMPDIR/manta_cram/workspace/pyflow.data/logs/pyflow_log.txt" ]; then
            echo "--- Manta CRAM log (last 30 lines) ---"
            tail -30 "$TMPDIR/manta_cram/workspace/pyflow.data/logs/pyflow_log.txt"
        fi
    else
        echo "configManta failed on CRAM — check output above"
    fi
else
    echo "Manta is not installed in this image."
    echo "To test, run this BAM through your Manta pipeline separately."
    echo ""
    echo "--- Simulating what Manta checks ---"
    echo "Proper pairs in BAM:  $(samtools view -f 0x2 -c "$TMPDIR/region.bam")"
    echo "MAPQ > 0 reads:      $(samtools view "$TMPDIR/region.bam" | head -n "$SAMPLE_N" | awk '$5 > 0' | wc -l | tr -d ' ') / $SAMPLE_N"
    echo "MAPQ distribution:"
    samtools view "$TMPDIR/region.bam" | head -n "$SAMPLE_N" \
        | awk '{mapq[$5]++} END {for(q in mapq) printf "  MAPQ=%s: %d reads\n", q, mapq[q]}' \
        | sort -t= -k2 -rn | head -10
fi

# ─────────────────────────────────────────────────────────────────────
echo ""
echo "=== Check 4: CRAM vs BAM readability comparison ==="
# ─────────────────────────────────────────────────────────────────────
echo "Reading 100 reads from CRAM with provided reference..."
CRAM_READS=$(samtools view --reference "$REF" "$CRAM" "$REGION" 2>/dev/null | head -100 | wc -l | tr -d ' ')
echo "  CRAM reads decoded: $CRAM_READS"

echo "Reading 100 reads from BAM (no reference needed)..."
BAM_READS=$(samtools view "$TMPDIR/region.bam" 2>/dev/null | head -100 | wc -l | tr -d ' ')
echo "  BAM reads decoded: $BAM_READS"

echo "Reading 100 reads from CRAM WITHOUT reference (simulates wrong/missing ref)..."
CRAM_NOREF=$(samtools view "$CRAM" "$REGION" 2>"$TMPDIR/noref_stderr" | head -100 | wc -l | tr -d ' ')
echo "  CRAM reads without ref: $CRAM_NOREF"
if [ -s "$TMPDIR/noref_stderr" ]; then
    echo "  stderr (first 5 lines):"
    head -5 "$TMPDIR/noref_stderr" | sed 's/^/    /'
fi

# ─────────────────────────────────────────────────────────────────────
echo ""
echo "=== Check 5: samtools version ==="
# ─────────────────────────────────────────────────────────────────────
samtools --version | head -3

# ─────────────────────────────────────────────────────────────────────
echo ""
echo "=== SUMMARY ==="
# ─────────────────────────────────────────────────────────────────────
echo "CRAM encoded with: Homo_sapiens_assembly38_masked.fasta (from @PG header)"
echo "Reference provided: $REF_BASENAME"
echo "CRAM reads with reference: $CRAM_READS"
echo "CRAM reads without reference: $CRAM_NOREF"
echo "BAM reads (self-contained): $BAM_READS"
echo ""
if [ "$CRAM_READS" -gt 0 ] && [ "$CRAM_NOREF" -eq 0 ]; then
    echo "KEY FINDING: CRAM requires the correct reference to decode."
    echo "If Manta is given a different reference (e.g. unmasked), it will see 0 reads."
    echo "Verify which reference your Manta pipeline uses."
elif [ "$CRAM_READS" -eq 0 ]; then
    echo "KEY FINDING: CRAM cannot be decoded even with the provided reference."
    echo "The reference may not match the encoding reference."
fi
echo ""
echo "NEXT STEPS:"
echo "1. Check which reference FASTA your Manta pipeline passes to Manta"
echo "2. If it differs from Homo_sapiens_assembly38_masked.fasta, that is the problem"
echo "3. Alternatively, convert to BAM first: samtools view -b --reference REF -o out.bam in.cram"

# Copy BAM to output path if requested
if [ -n "$OUTPUT_BAM" ]; then
    echo ""
    echo "=== Writing output BAM ==="
    cp "$TMPDIR/region.bam" "$OUTPUT_BAM" 2>/dev/null \
        || gsutil cp "$TMPDIR/region.bam" "$OUTPUT_BAM" 2>/dev/null \
        || echo "Failed to write output BAM to $OUTPUT_BAM"
    cp "$TMPDIR/region.bam.bai" "${OUTPUT_BAM}.bai" 2>/dev/null \
        || gsutil cp "$TMPDIR/region.bam.bai" "${OUTPUT_BAM}.bai" 2>/dev/null \
        || true
    echo "Wrote: $OUTPUT_BAM"
fi