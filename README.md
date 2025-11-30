# Bioinformatics – Primitive Local Alignment Search Tool (PLAST)

This repository contains my solution for **TP2 de bio-informatique (IFT3295, Université de Montréal)**.

The goal is to implement a simplified version of **BLASTN** called **PLAST** that:
- searches a query DNA sequence inside a database of sequences (FASTA),
- finds high-scoring local matches using a **seed-and-extend** strategy,
- extends High Scoring Pairs (HSPs),
- reports a **bitscore** and **E-value** for each significant hit.

Only one script is meant to be executed: `plast.py`.  
The rest of the repository contains the written report and support files for the assignment.

---

## How it works 

Given:
-i      query nucleotide sequence (string of A/C/G/T)
-db     FASTA file containing the sequence database
-E      greedy HSP extension cutoff (higher = longer extensions)
-ss     significance threshold on the E-value
-seed   seed pattern to use:
          - contiguous: '11111111111'
          - spaced: e.g. '111010010100110111' (PatternHunter-style)
          
the program:

1. **Builds a seed index** on the database  
   - Supports contiguous seeds (e.g. `'11111111111'`)  
   - And spaced seeds à la **PatternHunter** (e.g. `'111010010100110111'`)

2. **Finds seed hits** between the query and database sequences.

3. **Extends each hit** into a **High Scoring Pair (HSP)** using a greedy extension with cutoff `-E`.

4. **Scores each HSP** and computes:
   - raw **alignment score**,
   - **bitscore**,
   - **E-value**.

5. **Filters hits** by E-value threshold `-ss` and prints the best alignment(s) for each database sequence.

---

## Requirements

- Python 3.8+  

---

## How to run

From a terminal in the repository directory:

### Basic example

```bash
python plast.py -i CGTAGTCGGCTAACGCATACGCTTGATAAGCGTAAGAGCCC -db tRNAs.fasta -E 5 -ss 1e-3 -seed '11111111111'
