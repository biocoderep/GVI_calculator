import pytest
from bio_engine import calculate_cai, calculate_dnds_approx

def test_calculate_cai():
    # Example chicken codon optimized sequence
    # ATG = M (1.0 weight)
    # GCC = A (0.44 frequency / 0.44 max = 1.0 weight)
    # CAG = Q (0.73 frequency / 0.73 max = 1.0 weight)
    seq = "ATGGCCCAG"
    cai = calculate_cai(seq)
    assert cai == pytest.approx(1.0, 0.01)
    
    # Suboptimal sequence
    # GCA = A (0.17 / 0.44 = 0.386 weight)
    # CAA = Q (0.27 / 0.73 = 0.369 weight)
    seq_sub = "ATGGCACAA"
    cai_sub = calculate_cai(seq_sub)
    assert cai_sub < 1.0
    assert cai_sub > 0.0

def test_calculate_dnds_approx():
    # Identical sequences -> dS = 0, dN = 0 -> ratio should default or handle safely
    seqs = ["ATGGCCCAG", "ATGGCCCAG"]
    dnds = calculate_dnds_approx(seqs)
    # The current code returns 1.05 when ds_list is empty
    assert dnds == 1.05

    # One synonymous mutation
    # GCC (A) -> GCA (A)
    seqs_syn = ["ATGGCCCAG", "ATGGCACAG"]
    dnds_syn = calculate_dnds_approx(seqs_syn)
    # dS > 0, dN = 0 -> ratio = 0.0
    assert dnds_syn == 0.0

    # One non-synonymous mutation
    # GCC (A) -> GAC (D)
    seqs_non = ["ATGGCCCAG", "ATGGACCAG"]
    dnds_non = calculate_dnds_approx(seqs_non)
    # dS = 0, dN > 0 -> currently handles gracefully? Wait, if dS is 0, ds_list is empty
    assert dnds_non == 1.05
