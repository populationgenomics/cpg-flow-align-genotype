"""Derive candidate QC thresholds from every dataset's latest MultiQC report.

Run as two CPG Flow stages - see `qc_calibration_stages.py` and README.md in this
package. Every module here is pure: no Hail Batch, no stage machinery, so the analysis is
testable without either.
"""
