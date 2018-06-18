import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster

from pdb_io import parse_PDB, position_model
from q_value import calc_Q_value
from geometry import make_vectors_cone, calculate_angle, get_PCA_vector


def calc_helix_geometry(model, chain, start, end):
	"""
	Calculate Helix direction and position from chain, start and end residue id
	Returns: helix_vector, helix_center
	"""
	helix_CO_bonds = []
	helix_CAs = []
	if end - start + 1 < 5:
		return None, None
	else:
		# Save C=O bonds and CA atoms
		for i in range(start, end + 1):
			try:
				helix_CO_bonds.append(model[chain][i]['O'].coord - model[chain][i]['C'].coord)
				helix_CAs.append(model[chain][i]['CA'].coord)
			except KeyError:
				continue

	# Remove outliers via angle to sum vector
	sum_vector = sum(helix_CO_bonds)
	deviation_angles = [calculate_angle(sum_vector, bond, True) for bond in helix_CO_bonds]

	angle_mean = np.mean(deviation_angles)
	angle_std = np.std(deviation_angles)

	final_helix_CO_bonds = [helix_CO_bonds[i] for i in range(len(helix_CO_bonds))
							if abs(deviation_angles[i] - angle_mean) < angle_std]

	helix_vector = sum(final_helix_CO_bonds)
	helix_vector = helix_vector / np.linalg.norm(helix_vector)
	helix_center = sum(helix_CAs) / (1.0 * len(helix_CAs))

	return helix_vector, helix_center


def get_max_sliding_window(scores, window=15):
	"""
	Apply sliding window averaging over a list of scores
	Returns: Maximum score, start index of maximum window, list of window averages
	"""
	avgs = [sum(scores[i:(i + window)]) / window for i in range(len(scores) - window)]

	max_score = max(avgs)
	max_index = np.argmax(avgs)

	return max_score, max_index, avgs


def make_Q_value_plot(pre_slice_scores, pre_z_min, pre_z_max, pre_max_index, pre_max_score, vector,
					  best_slice_scores, best_z_min, best_z_max, best_max_index, best_max_score, best_vector,
					  window_size):
	"""
	Plots the Q value distribution of a model for a certain orientation before and after optimisation.
	Returns: figure object
	"""
	fig, ax = plt.subplots(figsize=(10, 7))
	ax.plot(range(pre_z_min, pre_z_max), pre_slice_scores, color='blue')
	ax.plot(range(best_z_min, best_z_max), best_slice_scores, color='green')
	ax.legend(["TM helix avg:[%1.3f, %1.3f, %1.3f] (%1.3f)" % (vector[0], vector[1], vector[2], pre_max_score),
			   "optimised:[%1.3f, %1.3f, %1.3f] (%1.3f)" % (best_vector[0], best_vector[1], best_vector[2], best_max_score)])
	ax.set_xlabel('Cross Section')
	ax.set_ylabel('Q-score')

	ax.plot([pre_z_min + pre_max_index, pre_z_min + pre_max_index + window_size], [pre_max_score, pre_max_score],
			color="#5555ff")
	ax.plot([best_z_min + best_max_index, best_z_min + best_max_index + window_size], [best_max_score, best_max_score],
			color="#66ff44")

	return fig


def optimise_vector(model, vector, window_size=30, n_opti=20):
	"""
	Optimises a membrane orientation vector for a given vector using the Q value objective function
	n_opti determines how many vectors should be tested within a 15° cone of the original vector
	Returns: Vector of highest score,
			 highest score,
			 z-coordinate of highest sliding window (middle of membrane)
			 figure depicting Q values before and after optimisation
	"""
	print("Optimising orientation... ", end='', flush=True)
	test_vectors = make_vectors_cone(vector, 15, n_opti)

	best_vector = vector
	best_slice_scores, best_z_min, best_z_max = calc_Q_value(model, vector)
	best_score, best_max_index, best_avgs = get_max_sliding_window(best_slice_scores, window_size)

	pre_slice_scores, pre_z_min, pre_z_max = best_slice_scores, best_z_min, best_z_max
	pre_max_index, pre_max_score = best_max_index, best_score

	for i, test_vector in enumerate(test_vectors):
		test_slice_scores, test_z_min, test_z_max = calc_Q_value(model, np.array(test_vector))
		test_score, test_max_index, test_avgs = get_max_sliding_window(test_slice_scores, window_size)

		if test_score > best_score:
			best_vector = test_vector
			best_score = test_score
			best_z_min = test_z_min
			best_z_max = test_z_max
			best_max_index = test_max_index
			best_slice_scores = test_slice_scores

	print("Complete!")
	fig = make_Q_value_plot(pre_slice_scores, pre_z_min, pre_z_max, pre_max_index, pre_max_score, vector,
							best_slice_scores, best_z_min, best_z_max, best_max_index, best_score, best_vector,
							window_size)

	z_opti = best_z_min + best_max_index + 0.5 * window_size
	return best_vector, best_score, z_opti, fig, pre_max_score


def process_results_file(filepath_to_csv, filepath_to_pdb, output_filepath, plot, n_opti):
	classi = pd.read_csv(filepath_to_csv)
	model = parse_PDB(filepath_to_pdb)
	model, init_pos = position_model(model)

	# Preprocess results-file
	vectors = []
	centers = []
	labels = []
	for i in classi.index:
		vector, center = calc_helix_geometry(model, classi.loc[i, "Chain"],
											 classi.loc[i, "Helix Start"], classi.loc[i, "Helix End"])
		vectors.append(vector)
		centers.append(center)
		labels.append(classi.loc[i, "Transmembrane Prediction"])

	helix_df = pd.DataFrame({'helix': classi.index, 'label': labels, 'vector': vectors, 'center': centers}).dropna()
	helix_df_tm = helix_df[np.array(helix_df['label']) == 1]


	print("Filtering false likely positive helices...")
	# Filter 1: by angle
	pca_vector1 = get_PCA_vector(np.array(helix_df_tm['vector']))
	angles = [(1 - np.abs(np.cos(calculate_angle(pca_vector1, vector)))) for vector in helix_df_tm['vector']]

	helix_df_tm = helix_df_tm.assign(angle_dist=angles, index=helix_df_tm.index)

	filtered1 = helix_df_tm[helix_df_tm['angle_dist'] > np.deg2rad(20)]
	if(len(filtered1.index) > 0):
		print("TM helices filtered by angle -> might be false positives")
		print(classi.loc[filtered1.index, ["Chain", "Helix Start", "Helix End"]])

	helix_df_tm = helix_df_tm[helix_df_tm['angle_dist'] <= np.deg2rad(20)]  # ODER FILTER 0.3 ~ 17.5°
	vectors_filter1 = [v.tolist() for v in helix_df_tm['vector']]
	centers_filter1 = [c.tolist() for c in helix_df_tm['center']]


	# Filter 2: by distance
	pca_vector2 = get_PCA_vector(np.array(helix_df_tm['vector']))
	projected_distances = [np.dot(pca_vector2, center) for center in helix_df_tm['center']]

	dist_matrix = pdist([[d, d] for d in projected_distances])
	Z = linkage(dist_matrix, 'single')
	clusters = fcluster(Z, 20, criterion='distance')
	biggest_cluster = max(clusters, key=clusters.tolist().count)

	filtered_indices = helix_df_tm.index[np.where(clusters != biggest_cluster)]
	filtered2 = helix_df_tm.loc[filtered_indices, :]
	if(len(filtered2.index) > 0):
		print("TM helices filtered by distance -> might be false positives")
		print(classi.loc[filtered2.index, ["Chain", "Helix Start", "Helix End"]])

	filter_indices = helix_df_tm.index[np.where(clusters == biggest_cluster)]
	helix_df_tm = helix_df_tm.loc[filter_indices, :]
	vectors_filter2 = [v.tolist() for v in helix_df_tm['vector']]
	centers_filter2 = [c.tolist() for c in helix_df_tm['center']]


	# Get final Orientation
	pca_vector3 = get_PCA_vector(np.array(helix_df_tm['vector']))
	center_final = sum(helix_df_tm['center']) / len(helix_df_tm['center'])
	print("Tentative orientation calculated.")
	best_vector, best_score, z_opti, fig, pre_max_score = optimise_vector(model, pca_vector3, window_size=30, n_opti=n_opti)

	print("Average TM helix vector: [%1.3f, %1.3f, %1.3f]"%(pca_vector3[0], pca_vector3[1], pca_vector3[2]) +
		  ", Score: %1.3f"%(pre_max_score))
	print("Optimised orientation: [%1.3f, %1.3f, %1.3f]"%(best_vector[0], best_vector[1], best_vector[2]) +
		  ", Score: %1.3f"%(best_score))

	if plot:
		fig.savefig(output_filepath.replace(".json", ".png"))

	results = {'init_pos': init_pos.tolist(), 'vectors_all': [v.tolist() for v in vectors],
			   'centers_all': [c.tolist() for c in centers], 'labels_all': labels,
			   'pca1': pca_vector1.tolist(),
			   'vectors_filter1': vectors_filter1, 'centers_filter1': centers_filter1,
			   'vectors_filter2': vectors_filter2, 'centers_filter2': centers_filter2,
			   'vector_final': pca_vector3.tolist(), 'center_final': center_final.tolist(),
			   'vector_opti': best_vector, 'z_opti': z_opti,
			   'score_before': pre_max_score, 'score_opti':best_score}

	with open(output_filepath, 'w') as file:
		json.dump(results, file)
	return results
