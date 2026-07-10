import os
from bio_engine import process_fasta_file
import json

filepath = "/mnt/e/Inlfuenza/avian_influenza/GVI_index_calculator/uploads/3c1c2ff4-fa10-493d-b17f-eb49c7d6c0a8_FMDV_2011_2015_aligned.fasta"

try:
    print(f"Running pipeline on {filepath}...")
    res = process_fasta_file(filepath)
    if res:
        print("\nPIPELINE FINISHED SUCCESSFULLY!")
        print(f"Detected Model: {res['Detected_Model']}")
        print(f"Evolutionary Rate: {res['Evolutionary_Rate']}")
        print(f"Recombination Rate: {res['Recombination_Rate']}")
        print(f"GVI Score: {res['GVI_Score']}%")
        print(f"Newick Tree Output snippet: {res['Newick_Tree'][:100]}")
    else:
        print("Pipeline returned None.")
except Exception as e:
    import traceback
    traceback.print_exc()
