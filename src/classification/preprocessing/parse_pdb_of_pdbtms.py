# Parse helices
import sys
import pandas as pd
import numpy as np
from prody import *
import random
import csv
import os
import urllib.request
print("Imports done")

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

helix_sequences = []
def parse_helices_pdb(pdb):
    atoms, header = parsePDB(pdb, header=True, )
    helix_ranges = header.get('helix_range')
    helix_indices = []
    for helix_range in helix_ranges:
        helix_chain = helix_range[1]
        helix_indices.append([helix_chain, helix_range[4],helix_range[5]])

    sequence = {}
    for atom in atoms:
        if atom.getChid() in sequence.keys():
            sequence.get(atom.getChid()).update({atom.getResnum() : d.get(atom.getResname())})
        else:
            sequence[atom.getChid()] = {atom.getResnum() : d.get(atom.getResname())}
    for helix_index in helix_indices:
        helix_chain = helix_index[0]
        helix_start_coords = []
        helix_end_coords = []
        helix_sequence = ""
        helix_start = helix_index[1]
        helix_end = helix_index[2]
        for key in range(helix_start,helix_end+1):
            helix_sequence = helix_sequence+sequence.get(helix_index[0]).get(key)
        helix_start_coords = list(atoms[helix_chain, helix_start]["CA"].getCoords())
        helix_end_coords = list(atoms[helix_chain, helix_end]["CA"].getCoords())
        helix = [header.get("identifier"), helix_chain, helix_start, helix_end, helix_start_coords, helix_end_coords, helix_end-helix_start+1, helix_sequence, "NaN"]
        helix_sequences.append(helix)
    return helix_sequences

# issues = []
# path="./pdbtm_pdb/"
# #path="./"
# for file in os.listdir(path):
#     if file.endswith(".gz"):
#         try:
#             parse_helices_pdb(path+file)
#         except:
#             print("ERROR:", file)
#             continue
#
#
# columns = ["PDB ID", "Chain", "Helix Start", "Helix End", "Helix Start CA", "Helix End CA", "Helix Length", "Helix Sequence", "Is Transmembrane"]
# df_helices = pd.DataFrame(helix_sequences, columns=columns)
# df_helices.to_csv("test_helices.csv", sep=',', encoding='utf-8', index=False)

