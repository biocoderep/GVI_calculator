import os
import subprocess
import logging
from Bio import SeqIO
from Bio.Seq import Seq

logger = logging.getLogger(__name__)

def parse_fasta(filepath):
    """Parses a FASTA file and validates it."""
    try:
        records = list(SeqIO.parse(filepath, "fasta"))
        if not records:
            logger.warning(f"No records found in {filepath}")
            return None
        return records
    except Exception as e:
        logger.error(f"Error parsing FASTA {filepath}: {e}")
        return None

def sanitize_for_phipack(fasta_path):
    """Sanitizes sequences for PhiPack by replacing unknown chars with '-'."""
    records = parse_fasta(fasta_path)
    if not records or len(records) < 4:
        return None 
    
    clean_records = []
    length = len(records[0].seq)
    valid_chars = set("ACTG-")
    
    for r in records:
        seq_str = str(r.seq).upper()
        cleaned_seq = "".join([c if c in valid_chars else "-" for c in seq_str])
        if len(cleaned_seq) == length:
            r.seq = Seq(cleaned_seq)
            clean_records.append(r)
            
    if len(clean_records) < 4:
        return None
        
    temp_path = fasta_path + ".clean.fasta"
    SeqIO.write(clean_records, temp_path, "fasta")
    return temp_path
