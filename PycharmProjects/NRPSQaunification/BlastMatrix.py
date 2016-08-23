import TwoSeqCompare

alt_dir = "Alternative Analysis"
ma_dir = "Matrices"

def create_blast_matrix(csv):
    for prime_seq in csv["non-specific sequence"]:
        for comp_seq in csv["non-specific sequence"]:
            TwoSeqCompare.alt_compare(prime_seq, comp_seq)
