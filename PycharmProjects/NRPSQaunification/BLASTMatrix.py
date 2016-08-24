import TwoSeqCompare
import os

alt_dir = "Alternative Analysis"
ma_dir = "Matrices"


def create_blast_matrix(csv):
    prime_list = []
    for prime_seq in csv:
        for comp_seq in reversed(csv):
                if not comp_seq in prime_list:
                    TwoSeqCompare.alt_compare(prime_seq, comp_seq)
        prime_list.append(prime_seq)



# creates a directory that does sit inside a non-major sub-directory
def create_dir(file, dir):
    if not os.path.exists(os.path.join(alt_dir, os.path.join(dir, file))):
        os.makedirs(os.path.join(alt_dir, os.path.join(dir, file)))


# creates a directory that does not sit inside a non-major sub-directory
def simple_dir(d):
    if not os.path.exists(os.path.join(alt_dir, d)):
        os.makedirs(os.path.join(alt_dir, d))
