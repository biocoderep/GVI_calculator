import sys
sys.path.insert(0, "/mnt/e/Inlfuenza/avian_influenza/GVI_index_calculator")
from scientific_engine.external_tools import extract_dates_from_fasta
d = extract_dates_from_fasta("/mnt/e/Inlfuenza/avian_influenza/GVI_index_calculator/uploads/3c1c2ff4-fa10-493d-b17f-eb49c7d6c0a8_FMDV_2011_2015_aligned.fasta")
print(d)
