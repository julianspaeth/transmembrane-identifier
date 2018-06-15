from Bio.PDB import *
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import numpy as np

import warnings


def parse_PDB(filepath):
    """
    Parses a PDB file from a given filepath using PDB Parser.
    Will check for Pilus & Virus proteins as well as DNA/RNA containing proteins.
    If a PDB file contains multiple models, the user is prompted to select one.
    :param filepath: Path to the PDB file
    :return: PDBParser model of the Protein
    """
    pdb_parser = PDBParser()

    # Actual parsing
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', PDBConstructionWarning)
        structure = pdb_parser.get_structure('protein', filepath)

    # Check for Pilus / Virus / DNA / RNA
    filter_words_name = ['pilus', 'pili', 'Pilus', 'Pili', 'PILUS', 'PILI', 'virus', 'Virus', 'VIRUS']
    filter_words_head = ['dna', 'DNA', 'rna', 'RNA']
    if any(word in structure.header['name'] for word in filter_words_name) or \
            any(word in structure.header['head'] for word in filter_words_head):
        print('Given structure is either a Pilus or Virus protein or contains nucleotides.')
        print('Will not be analysed.')
        return None

    # Select model from structure
    models = [model for model in structure]
    if len(models) > 1:
        model_number = input(str(len(models))+" models found. Please enter model id (1-index): ")
        model = models[int(model_number)-1]
    else:
        model = models[0]
    return model


def position_model(model):
    """
    Moves a PDBParser model to its center.
    :param model: PDBParser model
    :return: PDBParser model, with adjusted coordinates
    """
    total_vector = np.array([0, 0, 0])
    count = 0
    for chain in model:
        for residue in chain:
            for atom in residue:
                count += 1
                total_vector = total_vector + atom.coord
    previous_position = total_vector/count

    model.transform(rot=rotaxis2m(theta=0, vector=Vector(1, 0, 0)), tran=-previous_position)
    return model, previous_position

#model = parse_PDB('./structures/pdb1u9j.pdb')
#model = position_model(model)