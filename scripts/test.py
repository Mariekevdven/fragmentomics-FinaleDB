import pandas as pd
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

#--------------------------------------------------------------
#PARAMETERS
#--------------------------------------------------------------

data_dir = Path("/Users/mariekevandeven/stage_fragmentomics/data/")
MAPQ_cutoff = 55
max_len = 400

#---------------------------------------------------------
#SAMPLE PROCESSING FUNCTION
#---------------------------------------------------------

    # 1. histogram counts + overflow bin
def process_sample(path, MAPQ_cutoff=55, max_len=400):
    counts = np.zeros(max_len + 1, dtype=int)

    # 2. # fragments in raw data + after MAPQ filtering
    total_raw = 0
    total_filtered = 0

    # 3. open gzipped file met subprocess
    proc = subprocess.Popen(
        ["bgzip", "-cd", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    # 4. read output line by line, MAPQ filteren, lengte berekenen, histogram updaten
    for line in proc.stdout:
        total_raw += 1
        chr_n, start, end, mapQ, strand = line.strip().split("\t")

        mapQ = int(mapQ)

        #5. MAPQ filteren
        if mapQ < MAPQ_cutoff:
            continue    

        total_filtered += 1

        #6. lengte berekenen
        length = int(end) - int(start)

        #7. binning in histogram
        if length >= max_len:
            counts[max_len] += 1
        else:
            counts[length]  += 1

    #8. converteren naar pandas Series met overflow bin 
    raw_series = pd.Series(
        counts,
        index=list(range(max_len)) + [f"{max_len}+"] #f-string voor overflow bin (kan je makkelijk aanpassen)
    )

    # 9. normaliseren voor # gefilterde fragments
    norm_series = raw_series / total_filtered

    return raw_series, norm_series, total_filtered, total_raw

#---------------------------------------------------------
# END MOTIFS

#---------------------------------------------------------

chrom = chr_n

five_prime = str(genome[chrom][start:start+4]).upper()
three_prime = str(genome[chrom][end-4:end]).upper()

#---------------------------------------------------------
# MAIN PIPELINE
#---------------------------------------------------------

files = sorted(data_dir.glob("*.frag.tsv.bgz"))

print(f"Found {len(files)} samples")

#containers
raw_counts = {}
norm_counts = {}
qc_results = [] 

# loop
for file in files:
    sample_name = file.name.replace(".hg38.frag.tsv.bgz", "")

    print(f"Processing {sample_name}")

    raw, norm, total_filtered, total_raw = process_sample(file)

    #store results
    raw_counts[sample_name] = raw
    norm_counts[sample_name] = norm

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
qc_df = pd.DataFrame(qc_results)

# Bestanden opslaan
raw_df.to_csv("fragment_length_raw_counts.csv")
norm_df.to_csv("fragment_length_normalized_counts.csv")
qc_df.to_csv("fragment_length_qc_metrics.csv", index=False)

