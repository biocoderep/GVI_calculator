import time
import os
from scientific_engine.external_tools import calculate_beast_recombination, detect_best_evolutionary_model
fasta = "FMDV_2011_2015_aligned.fasta"
start = time.time()
model = detect_best_evolutionary_model(fasta)
print("Model:", model)
recomb, evo = calculate_beast_recombination(fasta, model)
print("Recomb:", recomb, "Evo:", evo)
print("Time taken:", time.time() - start)
