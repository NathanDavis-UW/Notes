import BLASTMatrix
import os

alt_root_dir = "Alternative NRPSRoot"

csv = []
for [dirpath, dirname, filename] in os.walk(alt_root_dir):
    csv.extend(filename)
BLASTMatrix.create_blast_matrix(csv)