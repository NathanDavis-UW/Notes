import TwoSeqCompare
import os
import MultiDimensionalScaler
import numpy as np
import pickle

alt_dir = "Alternative Analysis"
ma_dir = "Matrices"


def create_blast_matrix(csv):
    if os.path.exists("data_matrix_save.p") and os.path.exists("name_list_save.p"):
        name_list = pickle.load(open("name_list_save.p", "rb"))
        data_matrix = pickle.load(open("data_matrix_save.p", "rb"))
    else:
        prime_list = []
        name_list = []
        data_matrix = np.zeros((len(csv), len(csv)))
        i = 0
        for prime_seq in csv:
            name_list.append(prime_seq[:len(prime_seq) - 4])
            k = len(csv)-1
            for comp_seq in reversed(csv):
                if not comp_seq in prime_list:
                    print(data_matrix)
                    data_matrix.itemset((i, k), TwoSeqCompare.alt_compare(prime_seq, comp_seq))
                    data_matrix.itemset((k, i), data_matrix[i][k])
                k -= 1
            prime_list.append(prime_seq)
            i += 1
        pickle.dump(data_matrix, open("data_matrix_save.p", "wb"))
        pickle.dump(name_list, open("name_list_save.p", "wb"))
    MultiDimensionalScaler.create_mds(data_matrix, name_list)

# creates a directory that does sit inside a non-major sub-directory
def create_dir(file, dir):
    if not os.path.exists(os.path.join(alt_dir, os.path.join(dir, file))):
        os.makedirs(os.path.join(alt_dir, os.path.join(dir, file)))


# creates a directory that does not sit inside a non-major sub-directory
def simple_dir(d):
    if not os.path.exists(os.path.join(alt_dir, d)):
        os.makedirs(os.path.join(alt_dir, d))
