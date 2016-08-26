import TwoSeqCompare
import os
import MultiDimensionalScaler
import numpy as np

alt_dir = "Alternative Analysis"
ma_dir = "Matrices"


def create_blast_matrix(csv):
    prime_list = []
    name_list = []
    data_matrix = np.array([])
    for prime_seq in csv:

        name_list.append(prime_seq[:len(prime_seq) - 4])
        for comp_seq in csv:
            data_matrix = np.append(data_matrix, np.array([0]))
    np.reshape(data_matrix, (len(csv), len(csv)))
    i = 0
    for prime_seq in csv:
        k = len(csv)-1
        for comp_seq in reversed(csv):
                if not comp_seq in prime_list:
                    data_matrix[i][k] = TwoSeqCompare.alt_compare(prime_seq, comp_seq)
                    data_matrix[k][i] = data_matrix[i][k]
                k -= 1
        prime_list.append(prime_seq)
        i += 1
    MultiDimensionalScaler.create_mds(data_matrix, name_list)

# creates a directory that does sit inside a non-major sub-directory
def create_dir(file, dir):
    if not os.path.exists(os.path.join(alt_dir, os.path.join(dir, file))):
        os.makedirs(os.path.join(alt_dir, os.path.join(dir, file)))


# creates a directory that does not sit inside a non-major sub-directory
def simple_dir(d):
    if not os.path.exists(os.path.join(alt_dir, d)):
        os.makedirs(os.path.join(alt_dir, d))
