# Bio-BLASTn Alignment Filtering

This project provides a Python script that analyzes a circular DNA sequence to identify homologous regions within the sequence using BLASTn.

## Features

- Creates a local BLAST database from a single circular DNA FASTA sequence.
- Performs a BLASTn self-comparison to detect highly similar (>95%), non-overlapping regions.
- Filters alignments based on length and identity thresholds. Greedily takes the longest match.
- Outputs filtered homologous sequence alignments to a CSV file for downstream analysis.

## Requirements

- Python 3.x  
- biopython==1.83
- numpy==1.24.4

## Usage

1. Place your circular DNA sequence in a FASTA file (e.g., `NSCV1.fasta`) in the working directory.  
2. Run the script to generate a BLAST database, perform self-alignment, filter results, and save output as CSV.  
3. Check the output CSV file (`filtered_greedy_alignments_1.csv`) for the homologous alignments.
