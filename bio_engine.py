import os
import logging
from scientific_engine.parsers import parse_fasta
from scientific_engine.metrics import (
    calculate_cai, calculate_dnds_approx, calculate_pi_and_gd,
    calculate_mutation_burden, calculate_gc_content, calculate_gvi
)
from scientific_engine.external_tools import (
    calculate_beast_recombination, calculate_true_evolutionary_rate, calculate_bdsky_re, detect_best_evolutionary_model
)
import numpy as np
import job_manager

logger = logging.getLogger(__name__)

def process_fasta_file(filepath, custom_weights=None, job_id=None):
    """Main entry point for processing a FASTA file."""
    logger.info(f"Processing {filepath}")
    records = parse_fasta(filepath)
    # Aggressively sanitize sequence headers to prevent BEAST/IQ-TREE crashes
    seen_ids = set()
    invalid_chars = [',', '=', ':', '(', ')', '[', ']', ';', ' ', "'", '"']
    for r in records:
        for char in invalid_chars:
            r.id = r.id.replace(char, '_')
            r.name = r.name.replace(char, '_')
            r.description = r.description.replace(char, '_')
            
        original_id = r.id
        counter = 1
        while r.id in seen_ids:
            r.id = f"{original_id}_{counter}"
            counter += 1
        seen_ids.add(r.id)
        r.name = r.id
        r.description = r.id
        
    from Bio import SeqIO
    SeqIO.write(records, filepath, "fasta")
    
    sequences = [str(r.seq).upper() for r in records]
    
    if job_id: job_manager.update_job(job_id, progress=5, current_step="Step 1/9: Calculating Nucleotide Diversity (π)...")
    pi, gd = calculate_pi_and_gd(sequences)
    
    if job_id: job_manager.update_job(job_id, progress=15, current_step="Step 2/9: Calculating Mutation Burden (MB)...")
    mb = calculate_mutation_burden(sequences)
    
    if job_id: job_manager.update_job(job_id, progress=25, current_step="Step 3/9: Analyzing GC Content...")
    gc = calculate_gc_content(sequences)

    if job_id: job_manager.update_job(job_id, progress=35, current_step="Step 4/9: Calculating Codon Adaptation Index (CAI)...")
    mean_cai = np.mean([calculate_cai(seq) for seq in sequences])
    
    if job_id: job_manager.update_job(job_id, progress=45, current_step="Step 5/9: Estimating dN/dS Ratio...")
    mean_dnds = calculate_dnds_approx(sequences)

    # EXTERNAL WRAPPERS (TRUE ACCURACY)
    if job_id: job_manager.update_job(job_id, progress=55, current_step="Step 6/9: Running ModelFinder (Testing Substitution Models)...")
    detected_model = detect_best_evolutionary_model(filepath)

    if job_id: job_manager.update_job(job_id, progress=65, current_step="Step 7/9: Running BEAST MCMC (Recombination & Bayesian Clock)...")
    logger.debug(f"Running BEAST MCMC for {filepath}")
    true_recomb, beast_evo = calculate_beast_recombination(filepath, detected_model)
    
    if job_id: job_manager.update_job(job_id, progress=80, current_step="Step 8/9: Building IQ-TREE (LSD2 Clock & ML Tree)...")
    logger.debug(f"Running IQ-TREE for {filepath}")
    iqtree_evo, tree_newick = calculate_true_evolutionary_rate(filepath, detected_model)

    if job_id: job_manager.update_job(job_id, progress=90, current_step="Step 9/9: Running BDSKY MCMC (Reproduction Number)...")
    logger.debug(f"Running BDSKY MCMC for {filepath}")
    re_baseline = calculate_bdsky_re(filepath, detected_model)
    
    # We will use the highly rigorous BEAST Evolutionary Rate as the primary metric, 
    # but fallback to IQ-TREE's LSD2 if BEAST failed to sample it properly.
    final_evo_rate = beast_evo if beast_evo > 0.0000001 else iqtree_evo

    gvi_score = calculate_gvi(
        mu=final_evo_rate, 
        re=re_baseline, 
        pi=pi, 
        mb=mb, 
        dnds=mean_dnds, 
        gd=gd, 
        ri=true_recomb, 
        cai=mean_cai, 
        gc_content=gc,
        custom_weights=custom_weights
    )

    # Safely define year_str for output
    basename = os.path.basename(filepath)
    filename = basename.split("_", 1)[-1] if "_" in basename else basename
    year_str = filename.replace(".fasta", "").replace(".fa", "").replace("aligned_", "")

    results = {
        "Year": year_str,
        "Sequence_Count": len(records),
        "Detected_Model": detected_model,
        "Nucleotide_Diversity_Pi": round(pi, 5),
        "Genetic_Distance_GD": round(gd, 5),
        "Mutation_Burden_MB": round(mb, 2),
        "GC_Content": round(gc, 4),
        "Evolutionary_Rate": '{:.2e}'.format(final_evo_rate),
        "dN_dS_Ratio": round(mean_dnds, 4),
        "Codon_Adaptation_Index": round(mean_cai, 4),
        "Effective_Reproduction_Number_Re": round(re_baseline, 4),
        "Recombination_Rate": round(true_recomb, 4),
        "GVI_Score": round(gvi_score, 4),
        "Newick_Tree": tree_newick
    }
    
    logger.info(f"Finished processing {filepath}")
    return results
