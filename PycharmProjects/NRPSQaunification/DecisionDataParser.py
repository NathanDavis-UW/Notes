import os
import csvGenerator

root_dir = "NRPSRoot"


# parses genbank files to create a csv with abridged data for decision tree creation
def parse_data(root_analysis):
    files = []
    for filename in os.walk(root_dir):
        files.extend(filename[2])
    for f in files:
        if f[:len(f)-4] in root_analysis:
            similarities = []
            products = []
            length = 0
            organism = ""
            for line in f.readlines():
                if "note=" or "inference=" in line and "similar" or "similarity" in line :
                    similarities.append(line[line.index("=")+1:])
