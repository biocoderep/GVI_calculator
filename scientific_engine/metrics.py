import math
import itertools
import numpy as np

# --- BIOLOGICAL TABLES ---
STANDARD_GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

CHICKEN_CODON_USAGE = {
    'GCT': 0.28, 'GCC': 0.44, 'GCA': 0.17, 'GCG': 0.11, 'CGT': 0.11, 'CGC': 0.23, 'CGA': 0.11, 'CGG': 0.21, 'AGA': 0.19, 'AGG': 0.15,
    'AAT': 0.45, 'AAC': 0.55, 'GAT': 0.46, 'GAC': 0.54, 'TGT': 0.44, 'TGC': 0.56, 'CAA': 0.27, 'CAG': 0.73, 'GAA': 0.38, 'GAG': 0.62,
    'GGT': 0.21, 'GGC': 0.37, 'GGA': 0.24, 'GGG': 0.18, 'CAT': 0.45, 'CAC': 0.55, 'ATT': 0.41, 'ATC': 0.47, 'ATA': 0.12, 'TTA': 0.08,
    'TTG': 0.14, 'CTT': 0.14, 'CTC': 0.22, 'CTA': 0.08, 'CTG': 0.34, 'AAA': 0.41, 'AAG': 0.59, 'ATG': 1.00, 'TTT': 0.46, 'TTC': 0.54,
    'CCT': 0.28, 'CCC': 0.33, 'CCA': 0.26, 'CCG': 0.13, 'TCT': 0.18, 'TCC': 0.23, 'TCA': 0.14, 'TCG': 0.06, 'AGT': 0.15, 'AGC': 0.24,
    'ACT': 0.24, 'ACC': 0.43, 'ACA': 0.22, 'ACG': 0.11, 'TGG': 1.00, 'TAT': 0.44, 'TAC': 0.56, 'GTT': 0.21, 'GTC': 0.28, 'GTA': 0.14, 'GTG': 0.37
}

AMINO_ACID_MAX_WEIGHT = {}
for codon, freq in CHICKEN_CODON_USAGE.items():
    aa = STANDARD_GENETIC_CODE.get(codon)
    if aa and aa != '_':
        AMINO_ACID_MAX_WEIGHT[aa] = max(AMINO_ACID_MAX_WEIGHT.get(aa, 0), freq)

def calculate_cai(sequence):
    """Calculates Codon Adaptation Index based on Chicken codon usage."""
    valid_codons, cai_sum = 0, 0.0
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon in CHICKEN_CODON_USAGE:
            aa = STANDARD_GENETIC_CODE[codon]
            if aa != '_' and AMINO_ACID_MAX_WEIGHT[aa] > 0:
                w = CHICKEN_CODON_USAGE[codon] / AMINO_ACID_MAX_WEIGHT[aa]
                if w > 0:
                    cai_sum += math.log(w)
                    valid_codons += 1
    if valid_codons > 0: 
        return math.exp(cai_sum / valid_codons)
    return 0.0

def calculate_dnds_approx(sequences):
    """Calculates approximate dN/dS ratio for a set of sequences."""
    if len(sequences) < 2: return 0.0
    sample = sequences[:50] 
    ds_list, dn_list = [], []
    for s1, s2 in itertools.combinations(sample, 2):
        syn_mut, nonsyn_mut, valid_codons = 0, 0, 0
        for i in range(0, min(len(s1), len(s2)) - 2, 3):
            c1, c2 = s1[i:i+3], s2[i:i+3]
            if '-' in c1 or '-' in c2 or 'N' in c1 or 'N' in c2: continue
            if c1 != c2:
                aa1 = STANDARD_GENETIC_CODE.get(c1, 'X')
                aa2 = STANDARD_GENETIC_CODE.get(c2, 'X')
                if aa1 != 'X' and aa2 != 'X':
                    if aa1 == aa2: syn_mut += 1
                    else: nonsyn_mut += 1
            valid_codons += 1
        if valid_codons > 0:
            pS = syn_mut / (valid_codons * 0.25)
            pN = nonsyn_mut / (valid_codons * 0.75)
            if pS > 0 and pN > 0 and pS < 0.75 and pN < 0.75:
                dS = -0.75 * math.log(1 - (4/3 * pS))
                dN = -0.75 * math.log(1 - (4/3 * pN))
                if dS > 0:
                    ds_list.append(dS)
                    dn_list.append(dN)
    if ds_list:
        return sum(dn_list) / len(dn_list) / (sum(ds_list) / len(ds_list))
    return 1.05

def calculate_pi_and_gd(sequences):
    """Calculates Nucleotide Diversity (Pi) and Genetic Distance (GD)."""
    pi, gd = 0.0, 0.0
    if len(sequences) >= 2:
        diff_sum, comparisons = 0, 0
        for i in range(len(sequences)):
            for j in range(i+1, len(sequences)):
                diffs = sum(1 for a, b in zip(sequences[i], sequences[j]) if a != b and a != '-' and b != '-')
                valid = sum(1 for a, b in zip(sequences[i], sequences[j]) if a != '-' and b != '-')
                if valid > 0:
                    diff_sum += (diffs / valid)
                    comparisons += 1
        if comparisons > 0: 
            pi = gd = diff_sum / comparisons
    return pi, gd

def calculate_mutation_burden(sequences):
    """Calculates mutation burden relative to the first sequence."""
    if len(sequences) < 2: return 0.0
    ref = sequences[0]
    mb_list = []
    for seq in sequences[1:]:
        diffs = sum(1 for a, b in zip(ref, seq) if a != b and a != '-' and b != '-')
        valid = sum(1 for a, b in zip(ref, seq) if a != '-' and b != '-')
        if valid > 0: mb_list.append((diffs / valid) * 1000)
    return np.mean(mb_list) if mb_list else 0.0

def calculate_gc_content(sequences):
    """Calculates GC content average for the sequences."""
    gcs = []
    for seq in sequences:
        vseq = seq.replace("-", "").replace("N", "")
        if len(vseq) > 0: 
            gcs.append((vseq.count('G') + vseq.count('C')) / len(vseq))
    return np.mean(gcs) if gcs else 0.0

def calculate_gvi(mu, re, pi, mb, dnds, gd, ri, cai, gc_content, gc_reference=0.45, custom_weights=None):
    """
    Calculates the Genomic Virulence Index (GVI 2.0).
    Normalizes input parameters and computes the weighted sum.
    """
    # Normalization bounds (based on max thresholds from GVI 2.0 framework)
    try:
        mu_norm = min(1.0, max(0.0, float(mu) / 0.01))
    except (ValueError, TypeError):
        mu_norm = 0.0
        
    re_norm = min(1.0, max(0.0, re / 2.0))
    pi_norm = min(1.0, max(0.0, pi / 0.02))
    mb_norm = min(1.0, max(0.0, mb / 100.0))
    dnds_norm = min(1.0, max(0.0, dnds / 3.0))
    gd_norm = min(1.0, max(0.0, gd / 0.05))
    ri_norm = min(1.0, max(0.0, ri / 0.30))
    
    # GC deviation from a reference (assume 0.45 as baseline if not provided)
    gc_dev = abs(gc_content - gc_reference)
    gc_dev_norm = min(1.0, max(0.0, gc_dev / 0.03))
    
    # Combine CAI (0-1) and GC Deviation for w8
    cai_gc_norm = 0.5 * cai + 0.5 * gc_dev_norm
    
    # Default weights sum to 1.0
    w = [0.13, 0.30, 0.10, 0.20, 0.12, 0.10, 0.03, 0.02]
    
    if custom_weights:
        keys = ['w_mu', 'w_re', 'w_pi', 'w_mb', 'w_dnds', 'w_gd', 'w_ri', 'w_cai']
        for i, key in enumerate(keys):
            if key in custom_weights:
                w[i] = custom_weights[key]
                
    # Mathematically normalize array to sum to exactly 1.0 (100%)
    total_sum = sum(w)
    if total_sum > 0:
        w = [x / total_sum for x in w]
        # print(f"Normalized GVI Weights: {w} (Sum: {sum(w)})")
    
    gvi_score = (
        w[0] * mu_norm +
        w[1] * re_norm +
        w[2] * pi_norm +
        w[3] * mb_norm +
        w[4] * dnds_norm +
        w[5] * gd_norm +
        w[6] * ri_norm +
        w[7] * cai_gc_norm
    )
    
    return gvi_score
