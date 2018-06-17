import argparse

import pandas as pd
import numpy as np
from prody import parsePDB
from sklearn.externals import joblib


svm_path = 'svm.pkl'

def parse_arguments():
    """
    Parses the input arguments from the terminal
    """
    parser = argparse.ArgumentParser(description="Parse input for Transmembrane helix classifier.")
    parser.add_argument("-i", type=str, help="Path to pdb file", required=True)
    parser.add_argument("-o", type=str, help="Path to output CSV", required=True)
    return parser.parse_args()

def prepare_data(pdb_file):
    helix_df = parse_pdb_helices(pdb_file)
    print(helix_df)
    helix_df = helix_df.dropna()
    X = preprocess(helix_df)
    return X, helix_df

def preprocess(helix_df):
    print(helix_df)
    aas = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'H', 'I', 'G', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    def count_aas(helix):
        try:
            counts = []
            for aa in aas:
                counts.append(helix.count(aa))
        except:
            print("Error in preprocess the helix: counting amino acids was not possible")
        return counts

    def create_aa_frequency_features(df):
        X = []
        helices = list(df["Helix Sequence"])
        for helix in helices:
            X.append(count_aas(helix))
        return np.array(X)


    X = create_aa_frequency_features(helix_df)
    return X


def parse_pdb_helices(pdb_file):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    helix_sequences = []
    print(pdb_file)
    atoms, header = parsePDB("." + pdb_file, header=True)
    helix_ranges = header.get('helix_range')
    helix_indices = []
    for helix_range in helix_ranges:
        helix_chain = helix_range[1]
        helix_indices.append([helix_chain, helix_range[4],helix_range[5]])
    print(helix_indices)

    sequence = {}
    for atom in atoms:
        if atom.getChid() in sequence.keys():
            sequence.get(atom.getChid()).update({atom.getResnum() : d.get(atom.getResname())})
        else:
            sequence[atom.getChid()] = {atom.getResnum() : d.get(atom.getResname())}
    for helix_index in helix_indices:
        helix_sequence = ""
        helix_start = helix_index[1]
        helix_end = helix_index[2]
        helix_chain = helix_index[0]
        for key in range(helix_start,helix_end+1):
            try:
                helix_sequence = helix_sequence+sequence.get(helix_chain).get(key)
            except:
                continue
        try:
            helix_start_coords = list(atoms[helix_chain, helix_start]["CA"].getCoords())
            helix_end_coords = list(atoms[helix_chain, helix_end]["CA"].getCoords())
        except:
            helix_start_coords = list(atoms[helix_chain, helix_start]["CB"].getCoords())
            helix_end_coords = list(atoms[helix_chain, helix_end]["CB"].getCoords())
        helix_length = helix_end-helix_start+1
        helix = [header.get("identifier"), helix_chain, helix_start, helix_end, helix_start_coords, helix_end_coords, helix_length, helix_sequence, "NaN"]
        print(helix)
        if helix_length > 4:
            helix_sequences.append(helix)

    columns = ["PDB ID", "Chain", "Helix Start", "Helix End", "Helix Start CA", "Helix End CA", "Helix Length",
                   "Helix Sequence", "Transmembrane Prediction"]
    helix_sequences_df = pd.DataFrame(helix_sequences, columns=columns)
    return helix_sequences_df

def classify(X):
    svm = joblib.load(svm_path)
    y = svm.predict(X)
    return y

def create_output(path, df, y_pred):
    df['Transmembrane Prediction'] = y_pred
    df.to_csv(path, sep=',', encoding='utf-8', index=False)

if __name__ == "__main__":
    # Parse arguments
    print("Parse PDB file ...")
    args = parse_arguments()

    # Read pdb and preprocess for classification
    print("Preprocess the data ...")
    X, df = prepare_data(args.i)

    # Apply svm
    print("Classify helices ...")
    y = classify(X)

    # Create output CSV
    create_output(args.o, df, y)
    print("Finished!")


