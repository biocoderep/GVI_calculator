# GVI Calculator: Backend Algorithms and Biological Significance

This document outlines the specific algorithms and mathematical formulas utilized by the Genomic Vulnerability Index (GVI) Calculator backend. It details the technical calculations and their deep biological significance in the context of assessing viral threats such as Avian Influenza.

---

## 1. Codon Adaptation Index (CAI)

### The Formula
$$ CAI = \exp \left( \frac{1}{L} \sum_{i=1}^{L} \ln(w_i) \right) $$
*(Where $w_i$ is the relative adaptiveness of a codon compared to the most frequently used codon for that amino acid in the host).*

### Technical Logic
The backend compares the viral genetic sequence against the preferred codon usage table of a specific host (e.g., `CHICKEN_CODON_USAGE`). It calculates a score based on how perfectly the virus's codons match the host's preferences.

### Biological Importance
Viruses rely entirely on the host's cellular machinery (tRNA) to translate their RNA into proteins. A high CAI means the virus has genetically optimized itself to perfectly match the host. This results in **rapid, hyper-efficient viral replication** and a higher viral load, making the strain much more dangerous to that specific species.

---

## 2. Selection Pressure (dN/dS Ratio)

### The Algorithm
The engine uses a modified Jukes-Cantor calculation. It calculates the proportion of synonymous mutations ($pS$) and non-synonymous mutations ($pN$), applies a statistical correction, and divides the rates: 
$$ dS = -0.75 \ln(1 - \frac{4}{3}pS) $$
$$ dN = -0.75 \ln(1 - \frac{4}{3}pN) $$

### Technical Logic
It compares the rate of mutations that *change* the resulting protein ($dN$) against the rate of "silent" mutations that *do not change* the protein ($dS$). 

### Biological Importance
This is the ultimate signature of viral evolution. 
*   **$dN/dS < 1$ (Purifying Selection):** The virus is already highly successful; mutations are fatal to it, so it maintains its current form.
*   **$dN/dS > 1$ (Positive/Adaptive Selection):** The virus is actively mutating its proteins to survive. This is highly concerning, as it indicates the virus is actively adapting to evade host immune systems, vaccines, or antiviral drugs (Antigenic Drift).

---

## 3. Effective Reproductive Number ($R_e$) via BDSKY

### The Algorithm
**Birth-Death Skyline (BDSKY)** model running on Markov Chain Monte Carlo (MCMC) simulations via **BEAST 2**.

### Technical Logic
It mathematically infers the transmission rate of a disease *purely from genetic data*. By analyzing how rapidly the phylogenetic tree branches (the birth rate of new lineages vs. the death rate of extinct ones), the MCMC chain calculates $R_e$.

### Biological Importance
This reveals the real-time epidemic potential. If $R_e > 1$, the virus is spreading exponentially in the population. Calculating this directly from genomes (without relying on often-flawed hospital case data) gives a true, unbiased view of an ongoing outbreak.

---

## 4. Evolutionary Rate ($\mu$) & Molecular Clocks

### The Algorithm
Maximum Likelihood models (via IQ-TREE LSD2) or Bayesian MCMC relaxed clocks (via BEAST 2).

### Technical Logic
It calculates the "ticking speed" of the viral molecular clock—specifically, the number of nucleotide substitutions per site per year.

### Biological Importance
A high evolutionary rate means the virus is highly unstable and mutating rapidly. Pathogens with high evolutionary rates are dangerous because they can quickly generate novel variants that cross species barriers (zoonotic spillover) or rapidly bypass current vaccines.

---

## 5. Nucleotide Diversity ($\pi$)

### The Formula
Tajima's $\pi$ algorithm. 
$$ \pi = \frac{\sum d_{ij}}{n(n-1)/2} $$

### Technical Logic
It calculates the average number of nucleotide differences per site between any two random sequences in the uploaded dataset.

### Biological Importance
It measures the "genetic reservoir" of the viral population. High diversity means the virus exists as a massive swarm of slightly different variants. If a new environmental pressure is introduced (like a mass vaccination campaign), a highly diverse population is much more likely to contain a variant that is naturally immune to the vaccine.

---

## 6. Mutation Burden

### The Formula
Calculates the raw number of structural deviations from the oldest reference sequence, normalized per 1,000 base pairs.

### Technical Logic
It acts as a direct odometer of genetic drift, tracking exactly how far the currently circulating strains have drifted away from the original baseline strain.

### Biological Importance
High mutation burden strongly correlates with reduced vaccine efficacy. If the circulating strains have a massive mutation burden compared to the strain used to manufacture the vaccine, the vaccine antibodies will likely fail to recognize the new viruses.

---

## 7. GC Content & Deviation

### The Formula
$$ GC\% = \frac{G + C}{A + T + G + C} $$
*(The backend also calculates the deviation from a baseline stable reference, e.g., 0.45).*

### Technical Logic
Measures the structural proportion of Guanine and Cytosine bases in the viral RNA.

### Biological Importance
GC bonds are stronger than AT/AU bonds. Changes in GC content drastically alter the secondary 3D folding structure of the viral RNA. Viruses often alter their GC content to hide from the host's innate immune sensors (like Toll-Like Receptors) that specifically hunt for foreign RNA structures.

---

## Summary: The GVI Synthesis

By aggregating these parameters, the backend algorithm does not just look at *one* biological threat factor. It synthesizes **Transmission ($R_e$)**, **Adaptation (CAI & dN/dS)**, and **Instability ($\mu$ & $\pi$)** into a single mathematical index, providing a holistic, predictive score of how dangerous a specific viral lineage is becoming before a mass outbreak occurs.
