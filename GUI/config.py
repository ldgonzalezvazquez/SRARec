"""Application configuration for SRARec Desktop."""

APP_NAME = "SRARec Desktop"
APP_VERSION = "1.0.0"

# Replace with the real repository URL when available.
GITHUB_URL = "https://github.com/PLACEHOLDER/SRARec"

PAPER_STATUS = "Publication not yet available."

# Default analysis parameters (conservative)
DEFAULT_MIN_READS = 10
DEFAULT_MIN_GAMETE_READS = 5
DEFAULT_MIN_GAMETE_FREQ_PCT = 5.0

DEFAULT_PLATFORM = "ILLUMINA"
PLATFORMS = [
    ("Illumina", "ILLUMINA"),
    ("Oxford Nanopore", "OXFORD_NANOPORE"),
    ("PacBio", "PACBIO"),
]

REFERENCES = [
    ("SARS-CoV-2", "sars-cov-2"),
    ("HIV-1", "hiv-1"),
]

MIN_BASE_QUAL = 20

DEFAULT_OUTPUT_DIR_NAME = "SRARec_runs"
