import pandas as pd
import subprocess
import numpy as np
from pathlib import Path
from pyfaidx import Fasta
from itertools import product
import random
import sys
import argparse

#--------------------------------------------------------------
#ARGUMENT PARSER (for portability)
#--------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--reference-genome", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

#--------------------------------------------------------------
#PARAMETERS
#--------------------------------------------------------------

genome = Fasta(args.reference_genome)

MAPQ_cutoff = 55 #voor fragment quality filtering
valid_chroms = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
max_len = 500 #max fragment length voor histogram 
motif_sample_size = 100_000 #aantal fragments om te gebruiken voor motif counting (beperkt runtime)
random.seed(42)

#---------------------------------------------------------
# MOTIF SETUP
#---------------------------------------------------------

#alle 256 mogelijke 4-mers
motifs = [''.join(p) for p in product('ACGT', repeat=4)]

#reverse complement functie voor negatieve strand
def revcomp(seq):
    complement = str.maketrans('ACGT', 'TGCA')
    return seq.translate(complement)[::-1]

#---------------------------------------------------------
# subset selection en in functie gelijk fragment length feature extraction
#---------------------------------------------------------

    # 1. histogram counts + overflow bin
def process_sample(path, MAPQ_cutoff=55, max_len=500, k=100_000):
    counts = np.zeros(max_len + 1, dtype=int)

    # 2. drie dingen bouwen: total raw fragments, total filtered fragments, subset voor motif counting
    total_raw = 0
    mapq_filter = 0
    chrom_filter = 0
    total_filtered = 0
    subset = [] #voor motif counting

    # 3. open gzipped file met subprocess, lees line by line
    proc = subprocess.Popen(
        ["bgzip", "-cd", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    # 4. read output line by line, MAPQ filteren, lengte berekenen, histogram updaten en subset bijhouden
    for line in proc.stdout:
        total_raw += 1

        #6. MAPQ filteren 
        try: 
            chr_n, start, end, mapQ, strand = line.strip().split("\t")

            mapQ = int(mapQ)
            start = int(start)
            end = int(end)
        except ValueError:
            continue    

        if mapQ < MAPQ_cutoff:
            mapq_filter += 1
            continue

        if chr_n not in valid_chroms:
            chrom_filter += 1
            continue   

        total_filtered += 1 #onze echte data

        #fragment lengtes berekenen
        start = int(start)
        end = int(end)

        length = end - start

        #8. binning in histogram
        if length >= max_len:
            counts[max_len] += 1
        else:
            counts[length]  += 1

        # subset sampling voor motif counting (reservoir sampling)
        fragment = (chr_n, start, end, strand)

        #9. subset bijhouden, k = fragments voor motif counting
        if len(subset) < k:
            subset.append(fragment)
        else:
            j = random.randint(0, total_filtered - 1)
            if j < k:
                subset[j] = fragment

    # 10. converteren naar pandas Series met overflow bin 
    raw_series = pd.Series(
        counts,
        index=list(range(max_len)) + [f"{max_len}+"] #f-string voor overflow bin (kan je makkelijk aanpassen naar 500+ ofzo)
    )

    # 10. normaliseren voor # gefilterde fragments
    if total_filtered == 0:
        norm_series = pd.Series(0, index=raw_series.index)  # voorkom delen door nul
    else:
        norm_series = raw_series / total_filtered

    # 11. return output
    return raw_series, norm_series, total_filtered, total_raw, mapq_filter, chrom_filter, subset

        
#---------------------------------------------------------
# MOTIF COUNTING
#---------------------------------------------------------

# motif counting functie
def count_motifs(fragment_subset, genome):
    motif_counts = dict.fromkeys(motifs, 0)

    for chr_n, start, end, strand in fragment_subset:
        try:
            if start < 4 or end < 4:
                continue  # voor de zekerheid: skip fragments too close to chromosome ends 
                
            start_motif = genome[chr_n][start:start+4].seq.upper() #eerste 4 bases
            end_motif = revcomp(genome[chr_n][end-4:end].seq.upper()) #laatste 4 bases, reverse complement omdat we altijd in 5'->3' richting willen kijken

            #counting 
            if start_motif in motif_counts:
                motif_counts[start_motif] += 1

            if end_motif in motif_counts:
                motif_counts[end_motif] += 1

        except KeyError:
            continue 
    
    #normalisatie voor totaal aantal getelde motifs, als je niet 100K motifs over houdt na MapQ filtering
    total_motifs = sum(motif_counts.values())

    if total_motifs > 0:
        for motif in motif_counts:
            motif_counts[motif] /= total_motifs #output is distributie van motifs in de subset

    return motif_counts, total_motifs


#---------------------------------------------------------
# SLURM SCRIPT
#---------------------------------------------------------

# Open file from argparser
file = Path(args.input)
sample_name = file.name.replace(".hg38.frag.tsv.bgz", "") #hier dus alleen de sample naam

print(f"Processing sample {sample_name}")

#--------------------------------------------------------
# RUN SAMPLE PIPELINE
#--------------------------------------------------------

raw, norm, total_filtered, total_raw, mapq_filter, chrom_filter, subset = process_sample(
    file,
    MAPQ_cutoff=MAPQ_cutoff,
    max_len=max_len,
    k=motif_sample_size
    )

motif_counts, total_motifs = count_motifs(subset, genome)

#QC table met total raw, total filtered, ratio
qc = pd.DataFrame([{
    "sample": sample_name,
    "raw_fragments": total_raw,
    "mapq_filtered": mapq_filter, # How many fragments outside mapq filter?
    "chrom_filtered": chrom_filter, # How many fragments not on valid_chroms?
    "filtered_fragments": total_filtered,
    "lost_fragments": total_raw - total_filtered, #hoeveelheid weggefilterde fragments
    "total_motifs": total_motifs, # totaal aantal motifs (om motif fractions terug te rekenen naar total counts)
    "ratio": (
        total_filtered / total_raw) 
        if total_raw > 0 
        else np.nan # voorkom delen door nul
}])

#---------------------------------------------------------
# SAVE OUTPUT Per sample
#---------------------------------------------------------

outdir = Path(args.outdir)
outdir.mkdir(exist_ok=True)

raw.to_csv(outdir / f"{sample_name}_raw_length_distribution.csv")
norm.to_csv(outdir / f"{sample_name}_norm_length_distribution.csv")

pd.Series(motif_counts).to_csv(outdir / f"{sample_name}_motifs.csv")
qc.to_csv(outdir / f"{sample_name}_QC.csv", index=False)
