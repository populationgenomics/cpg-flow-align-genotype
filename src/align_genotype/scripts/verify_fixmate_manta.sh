#!/usr/bin/env bash
set -eu

# Verify that missing MQ/MC tags are the cause of Manta's
# "0 high-confidence read pairs" failure.
#
# Extracts a region, runs samtools fixmate, and compares
# MQ/MC tag counts before vs after.
#
# Usage:
#   ./verify_fixmate_manta.sh <sample.cram> <reference.fasta>

CRAM="${1:?Usage: $0 <sample.cram> <reference.fasta>}"
REF="${2:?Usage: $0 <sample.cram> <reference.fasta>}"
REGION="chr1:10000000-50000000"
SAMPLE_N=10000
CPU=3

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== Step 1: Extract region from CRAM ==="
samtools view -b --reference "$REF" -o "$TMPDIR/region.bam" "$CRAM" "$REGION"
samtools index "$TMPDIR/region.bam"
TOTAL=$(samtools view -c "$TMPDIR/region.bam")
echo "Extracted reads: $TOTAL"

echo ""
echo "=== Step 2: Check tags BEFORE fixmate ==="
ORIG_MQ=$(samtools view "$TMPDIR/region.bam" | head -n "$SAMPLE_N" | grep -c 'MQ:i:' || true)
ORIG_MC=$(samtools view "$TMPDIR/region.bam" | head -n "$SAMPLE_N" | grep -c 'MC:Z:' || true)
echo "MQ tags: $ORIG_MQ / $SAMPLE_N"
echo "MC tags: $ORIG_MC / $SAMPLE_N"

echo ""
echo "=== Step 3: Run samtools fixmate ==="
samtools sort -n -@ "$CPU" -o "$TMPDIR/namesorted.bam" "$TMPDIR/region.bam"
samtools fixmate -m "$TMPDIR/namesorted.bam" "$TMPDIR/fixed_namesorted.bam"
samtools sort -@ "$CPU" -o "$TMPDIR/fixed.bam" "$TMPDIR/fixed_namesorted.bam"
samtools index "$TMPDIR/fixed.bam"

echo ""
echo "=== Step 4: Check tags AFTER fixmate ==="
FIXED_MQ=$(samtools view "$TMPDIR/fixed.bam" | head -n "$SAMPLE_N" | grep -c 'MQ:i:' || true)
FIXED_MC=$(samtools view "$TMPDIR/fixed.bam" | head -n "$SAMPLE_N" | grep -c 'MC:Z:' || true)
echo "MQ tags: $FIXED_MQ / $SAMPLE_N"
echo "MC tags: $FIXED_MC / $SAMPLE_N"

echo ""
echo "=== Step 5: Sample reads comparison ==="
echo "--- Original (first proper pair) ---"
samtools view -f 0x2 "$TMPDIR/region.bam" | head -1 | tr '\t' '\n' | grep -E '^(MQ|MC|ms):' || echo "  (no MQ/MC/ms tags)"
echo "--- Fixed (first proper pair) ---"
samtools view -f 0x2 "$TMPDIR/fixed.bam" | head -1 | tr '\t' '\n' | grep -E '^(MQ|MC|ms):' || echo "  (no MQ/MC/ms tags)"

echo ""
echo "=== CONCLUSION ==="
echo "Original MQ: $ORIG_MQ / $SAMPLE_N"
echo "Fixed MQ:    $FIXED_MQ / $SAMPLE_N"
echo "Original MC: $ORIG_MC / $SAMPLE_N"
echo "Fixed MC:    $FIXED_MC / $SAMPLE_N"
if [ "$ORIG_MQ" -eq 0 ] && [ "$FIXED_MQ" -gt 0 ]; then
    echo ""
    echo "CONFIRMED: samtools fixmate adds MQ/MC tags to this CRAM."
    echo "Missing MQ/MC tags are the likely cause of Manta failure."
    echo "Fix: run 'samtools fixmate -m' on the full CRAM before Manta."
elif [ "$ORIG_MQ" -gt 0 ]; then
    echo ""
    echo "MQ tags were already present — MQ/MC is NOT the issue."
else
    echo ""
    echo "fixmate did not add MQ tags — investigate further."
fi