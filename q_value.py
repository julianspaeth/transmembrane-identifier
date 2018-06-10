from Bio.PDB import *
import freesasa
import numpy as np
import matplotlib.pyplot as plt

from pdb_io import *
from geometry import get_rotation


def get_atom_info():
	"""
	Returns a dictionary of atomic radii, hydrophobic and hydrophilic residues
	needed for the calculation
	:return: atom_info - Dictionary containing: atom_radii -> Dict: atom_radii['element-symbol'] = radius
												hydrophobic -> List of 3-letter hydrophobic residues
												hydrophilic -> List of 3-letter hydrophilic residues
	"""
	atom_radii = {'C': 1.80, 'N': 1.60, 'O': 1.40, 'S': 1.85}
	hydrophobic = ['PHE', 'GLY', 'ILE', 'LEU', 'MET', 'VAL', 'TRP', 'TYR']
	hydrophilic = ['ALA', 'CYS', 'ASP', 'GLU', 'ASN', 'GLN', 'HIS', 'LYS', 'PRO', 'ARG', 'SER', 'THR']
	atom_info = {'atom_radii': atom_radii, 'hydrophobic': hydrophobic, 'hydrophilic': hydrophilic}
	return atom_info


def get_z_bounds(model):
	"""
	Returns integer z-coordinate boundaries for a given PDBParser model
	:param model: PDBParser model
	:return: z_min - Integer - z-position of the lowest occurring atom
			 z_max - Integer - z-position of the highest occurring atom
	"""
	max_z = -100000
	min_z = 100000
	for chain in model:
		for residue in chain:
			for atom in residue:
				coordinates = atom.coord
				if coordinates[2] > max_z:
					max_z = coordinates[2]
				if coordinates[2] < min_z:
					min_z = coordinates[2]
	return int(np.floor(min_z)), int(np.ceil(max_z))


def make_slices(model, z_min, z_max):
	"""
	Splits a model into a dictionary of atom lists - separated by z-coordinates
	:param model: PDBParser model
	:param z_min: Integer - Lowest possible z-coordinate of an Atom, rounded down
	:param z_max: Integer - Highest possible z-coordinate of an Atom, rounded up
	:return: slices - Dictionary of atom lists, aggregated by z-coordinates
			 slices[10] = List of Atoms with z in 10-11 A
	"""
	slices = {}
	for i in range(z_min, z_max):
		slices[i] = []

	for chain in model:
		for residue in chain:
			for atom in residue:
				if atom.element in ['FE','H', 'MG']:
					continue
				slices[int(np.floor(atom.coord[2]))].append(atom)
	return slices


# HYDROPHOBIC FACTOR
def rectangle_grid(slice_list, z):
	"""
	Calculates a list of points that form a rectangle which encloses
	all atoms of the given slice at the height (z-coordinate)
	:param slice_list: List of Atoms - Slice to make rectangle for
	:param z: float - Height of the rectangle
	:return: rectangle_points: List of points
	"""
	max_x = -1000
	min_x = 1000
	max_y = -1000
	min_y = 1000
	for atom in slice_list:
		coordinates = atom.coord
		if coordinates[0] > max_x:
			max_x = coordinates[0]
		if coordinates[0] < min_x:
			min_x = coordinates[0]
		if coordinates[1] > max_y:
			max_y = coordinates[1]
		if coordinates[1] < min_y:
			min_y = coordinates[1]

	rectangle_points = []
	offset = 1

	for x in np.linspace(min_x - offset, max_x + offset, max([int(len(slice_list) / 3), 2])):
		rectangle_points.append(np.array([x, min_y - offset, z]))
		rectangle_points.append(np.array([x, max_y + offset, z]))
	for y in np.linspace(min_y - offset, max_y + offset, max([int(len(slice_list) / 3), 2])):
		rectangle_points.append(np.array([min_x - offset, y, z]))
		rectangle_points.append(np.array([max_x + offset, y, z]))
	return rectangle_points


def get_outer_atoms(slice_list, rectangle_points):
	outer_atoms = set()
	outer_atom_ids = []

	for point in rectangle_points:
		current_dists = []
		for atom in slice_list:
			current_dists.append(np.linalg.norm(point - atom.coord))

		outer_atoms.add(slice_list[np.argmin(current_dists)])

	for atom in list(outer_atoms):
		outer_atom_ids.append(atom.serial_number)

	return outer_atoms, outer_atom_ids


def calc_hydrophobic_factor(slice_list, slice_height):
	"""
	Calculates the hydrophobic factor of a given slice with a specified height
	:param slice_list: List of Atoms
	:param slice_height: Integer - Lower bound height (z-coordinate) of the current slice
	:return: hydrophobic factor - float
	"""
	slice_coords = []
	slice_radii = []
	atom_info = get_atom_info()

	for atom in slice_list:
		slice_radii.append(atom_info['atom_radii'][atom.element])
		slice_coords.extend(atom.coord)

	rectangle_points = rectangle_grid(slice_list, slice_height + 0.5)
	outer_atoms, outer_atom_ids = get_outer_atoms(slice_list, rectangle_points)

	sasa_result = freesasa.calcCoord(coord=slice_coords, radii=slice_radii)
	hydrophobic_sasa = 0
	total_sasa = 0

	for i, atom in enumerate(slice_list):
		if atom.serial_number in outer_atom_ids:
			if (atom.parent.resname in atom_info['hydrophobic']) or (atom.parent.resname in atom_info['hydrophilic']):
				total_sasa += sasa_result.atomArea(i)
				if atom.parent.resname in atom_info['hydrophobic']:
					hydrophobic_sasa += sasa_result.atomArea(i)
	if total_sasa == 0:
		return 0.00001
	else:
		return hydrophobic_sasa / total_sasa


# STRUCTURE FACTOR
def make_residue_dictionary(model):
	"""
	Makes a dictionary of CA-Atoms from a given model
	:param model: PDBParser model
	:return: residue_dict - Nested dictionary containing CA-atoms:
			 residue_dict[chain][resi] = <CA-Atom>
	"""
	residue_dict = {}

	for chain in model:
		residue_dict[chain.id] = {}
		for residue in chain:
			if residue.id[0] == ' ':
				residue_dict[chain.id][residue.id[1]] = residue['CA']
	return residue_dict


def get_unique_residues_from_slice(slice_list):
	"""
	Returns a list of unique residue IDs from a given slice list
	:param slice_list: List of Atoms
	:return: List of unique residue IDs
	"""
	ids = set()
	for atom in slice_list:
		if atom.parent.id[0] == ' ':
			ids.add((atom.parent.parent.id, atom.parent.id[1]))
	return list(ids)


def calc_structure_factor(slice_list, residue_dict, vector):
	"""
	Calculates the hydrophobic factor of a given slice with a specified height
	:param slice_list: List of Atoms
	:param residue_dict: Nested dictionary of CA-Atoms (make_residue_dictionary)
	:param vector: Normal vector of currently tested membrane plane
	:return: structure factor - float
	"""
	residue_ids = get_unique_residues_from_slice(slice_list)
	if len(residue_ids) == 0:
		return 0.00001
	straight_count = 0
	turn_count = 0
	endchain_count = 0

	for chain, resi in residue_ids:
		# Residue i-3
		try:
			prev_res = residue_dict[chain][resi - 3].coord
		except KeyError:
			endchain_count += 1
			continue

		# Residue i
		curr_res = residue_dict[chain][resi].coord

		# Residue i+3
		try:
			next_res = residue_dict[chain][resi + 3].coord
		except KeyError:
			endchain_count += 1
			continue

		prev_value = np.dot(vector, prev_res)
		curr_value = np.dot(vector, curr_res)
		next_value = np.dot(vector, next_res)

		if (prev_value < curr_value and curr_value < next_value) or (
				prev_value > curr_value and curr_value > next_value):
			straight_count += 1
		else:
			turn_count += 1

	straightness_factor = straight_count / len(residue_ids)
	turn_factor = 1 - (turn_count / len(residue_ids))
	endchain_factor = 1 - (endchain_count / len(residue_ids))

	return straightness_factor * turn_factor * endchain_factor


# Q value calculation
def calc_Q_value(model, vector):
	"""
	Calculates the Q-values of a given model for a given membrane plane normal vector
	:param model: PDBParser model
	:param vector: List of 3 floats - 3D vector of the tested normal vector
	:return: slice_scores - List of floats - Q-values for each slice along the vector
			 z_min - Integer - Lowest possible z-coordinate of the model along the vector, rounded down
			 z_max - Integer - Highest possible z-coordinate of the model along the vector, rounded up
	"""
	model_copy = model.copy()

	# Rotate model so that it aligns with [0,0,1]
	rot_axis, angle = get_rotation(vector)
	model_copy.transform(rot=rotaxis2m(theta=angle, vector=Vector(rot_axis)), tran=0)

	# get z-Coordinate bounds
	z_min, z_max = get_z_bounds(model_copy)
	# generate slices
	slices = make_slices(model_copy, z_min, z_max)
	# make dictionary of CA-atoms
	residue_dict = make_residue_dictionary(model_copy)

	# calculate Q-scores for each slice
	slice_scores = []
	for slice_height in slices.keys():
		slice_list = slices[slice_height]

		if len(slice_list) == 0:
			slice_scores.append(0)
			continue
		hydrophobic_factor = calc_hydrophobic_factor(slice_list, slice_height)
		structure_factor = calc_structure_factor(slice_list, residue_dict, np.array(vector))
		slice_scores.append(hydrophobic_factor * structure_factor)
		#print("Slice: %s\tAtoms: %s\tHydro_factor: %1.3f\tStruc_factor: %1.3f"%
			  #(slice_height, len(slice_list), hydrophobic_factor, structure_factor))

	return slice_scores, z_min, z_max






#model1 = parse_PDB('./structures/pdb1u9j.pdb')

#model1 = position_model(model1)
#z_min, z_max = get_z_bounds(model1)



