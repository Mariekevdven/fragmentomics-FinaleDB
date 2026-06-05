import pandas as pd
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from pyfaidx import Fasta
from itertools import product

#--------------------------------------------------------------
#PARAMETERS
#--------------------------------------------------------------

data_dir = Path("/Users/mariekevandeven/stage_fragmentomics/data/")
genome = Fasta("/Users/mariekevandeven/stage_fragmentomics/data/hg38.fa")
MAPQ_cutoff = 55
max_len = 400

#---------------------------------------------------------
# MOTIF SETUP
#---------------------------------------------------------

motifs = [''.join(p) for p in product('ACGT', repeat=4)]

def revcomp(seq):
    complement = str.maketrans('ACGT', 'TGCA')
    return seq.translate(complement)[::-1]

#---------------------------------------------------------
#SAMPLE PROCESSING FUNCTION
#---------------------------------------------------------

    # 1. histogram counts + overflow bin
def process_sample(path, MAPQ_cutoff=55, max_len=400):
    counts = np.zeros(max_len + 1, dtype=int)

    # 2. # fragments in raw data + after MAPQ filtering
    total_raw = 0
    total_filtered = 0

    # 3. motif counts initialiseren
    motif_counts_5 = dict.fromkeys(motifs, 0)
    motif_counts_3 = dict.fromkeys(motifs, 0)

    # 4. open gzipped file met subprocess
    proc = subprocess.Popen(
        ["bgzip", "-cd", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    # 5. read output line by line, MAPQ filteren, lengte berekenen, histogram updaten
    for line in proc.stdout:
        total_raw += 1
        chr_n, start, end, mapQ, strand = line.strip().split("\t")

        mapQ = int(mapQ)

        #6. MAPQ filteren
        if mapQ < MAPQ_cutoff:
            continue    

        total_filtered += 1

        start = int(start)
        end = int(end)

        #7. lengte berekenen
        length = int(end) - int(start)

        #8. binning in histogram
        if length >= max_len:
            counts[max_len] += 1
        else:
            counts[length]  += 1
        #----------------------------------------
        # End Motifs
        #----------------------------------------

        #9. chrom correctie (chr weghalen)
        try :
            chrom = chr_n

            #9.a correct for strand 
            if strand == "+":
                five_prime = str(genome[chrom][start:start+4]).upper()
                three_prime = str(genome[chrom][end-4:end]).upper()
        
            elif strand == "-":
                five_prime = revcomp(str(genome[chrom][start:start+4]).upper())
                three_prime = revcomp(str(genome[chrom][end-4:end]).upper())
            else:
                continue

            #9.b count motifs
            if five_prime in motif_counts_5:
                motif_counts_5[five_prime] += 1
            if three_prime in motif_counts_3:
                motif_counts_3[three_prime] += 1

        #9.c skip when chrom not found in genome 
        except KeyError:
            continue

    # 10. converteren naar pandas Series met overflow bin 
    raw_series = pd.Series(
        counts,
        index=list(range(max_len)) + [f"{max_len}+"] #f-string voor overflow bin (kan je makkelijk aanpassen)
    )

    # 10. normaliseren voor # gefilterde fragments
    norm_series = raw_series / total_filtered

    # 11. return output
    return raw_series, norm_series, total_filtered, total_raw, motif_counts_5, motif_counts_3

#---------------------------------------------------------
# MAIN PIPELINE
#---------------------------------------------------------

files = sorted(data_dir.glob("*.frag.tsv.bgz"))

print(f"Found {len(files)} samples")

#containers
raw_counts = {}
norm_counts = {}
qc_results = [] 

motif5_matrix = {}
motif3_matrix = {}

# loop
for file in files[:1]: #[:1] voor test run, haal dit weg voor volledige run
    sample_name = file.name.replace(".hg38.frag.tsv.bgz", "")

    print(f"Processing {sample_name}")

    raw, norm, total_filtered, total_raw, motif_counts_5, motif_counts_3 = process_sample(file)

    #store length counts
    raw_counts[sample_name] = raw
    norm_counts[sample_name] = norm

    #store motif counts
    motif5_matrix[sample_name] = motif_counts_5
    motif3_matrix[sample_name] = motif_counts_3

    #QC table met total raw, total filtered, ratio
    qc_results.append({
        "sample": sample_name,
        "raw_fragments": total_raw,
        "filtered_fragments": total_filtered,
        "ratio": total_filtered / total_raw
    })
# dataframes maken
raw_df = pd.DataFrame(raw_counts).T
norm_df = pd.DataFrame(norm_counts).T

motifs_5 = pd.DataFrame(motif5_matrix).T
motifs_3 = pd.DataFrame(motif3_matrix).T

qc_df = pd.DataFrame(qc_results)

final_df = pd.concat([norm_df, motifs_5, motifs_3], axis=1)

# Bestanden opslaan
raw_df.to_csv("fragment_length_raw_counts.csv")
norm_df.to_csv("fragment_length_normalized_counts.csv")
qc_df.to_csv("fragment_length_qc_metrics.csv", index=False)
motifs_5.to_csv("fragment_length_motifs_5.csv")
motifs_3.to_csv("fragment_length_motifs_3.csv")

final_df.to_csv("fragment_length_final_features.csv")
