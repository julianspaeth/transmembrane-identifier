import numpy as np
import math

from sklearn.decomposition import PCA


def get_rotation(vector):
	"""
	Computes the rotation axis and angle for a given 3D vector
	such that after the rotation, the vector will align with [0,0,1].
	:param vector: List of 3 numbers
	:return: rotation_axis: numpy array of the rotational axis
			 angle: Rotation angle in rad
	"""
	if np.linalg.norm(vector - np.array([0, 0, 1])) == 0:
		return np.array([1, 0, 0]), 0
	if np.linalg.norm(vector - np.array([0, 0, -1])) == 0:
		return np.array([1, 0, 0]), np.pi

	vector2 = np.array(vector)
	rotation_axis = np.cross(np.array([0, 0, 1]), vector2)

	angle = math.asin(np.linalg.norm(rotation_axis) / np.linalg.norm(vector2))
	return rotation_axis, angle


def calculate_angle(a, b, deg=False):
	"""
	Calculates the angle between two 3D vectors.
	:param a: Numpy array with 3 coordinates
	:param b: Numpy array with 3 coordinates
	:param deg: boolean - True -> result in Degrees, False -> result in rad
	:return: Number - angle in rad (or Degrees if deg=True)
	"""
	if a[0] == b[0] and a[1] == b[1] and a[2] == b[2]:
		return 0

	angle = math.acos(np.dot(a, b) / (np.linalg.norm(b) * np.linalg.norm(a)))
	if not deg:
		return angle
	else:
		return np.degrees(angle)


def make_vectors_sphere(n=100):
	"""
	Generates an array of evenly distributed 3D unit vectors using a Fibonacci sphere.
	Adapted from
		https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere/26127012#26127012
	:param n: Number of vectors to be generated
	:return: List of generated vectors
	"""
	vector_list = []
	y_offset = 2.0 / n
	step_size = math.pi * (3.0 - math.sqrt(5.0))

	for i in range(n):
		y = ((i * y_offset) - 1) + (y_offset / 2)
		r = math.sqrt(1 - y*y)

		phi = ((i + 1.0) % n) * step_size

		x = math.cos(phi) * r
		z = math.sin(phi) * r

		vector_list.append([x, y, z])

	return vector_list


def interpolate(b, xs, ys):
	"""
	Helper function for make_vectors_cone()
	"""
	for i, x in enumerate(xs):
		if b == x:
			return ys[i]
		if b < x:
			return ys[i-1] + (b-ys[i-1]) * (ys[i]-ys[i-1])/(xs[i]-xs[i-1])


def make_vectors_cone(vector, angle=180, desired_number=100):
	"""
	Generates an array of evenly distributed 3D unit vectors along a
	predefined vector with a maximal deviation of angle.
	:param vector: 3D numpy array
	:param angle: Number - Angle cutoff in Degrees
	:param desired_number: Number - Roughly how many vectors should be generated
	:return: List of generated vectors
	"""
	xs = [0, 1, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180]
	ys = [0, 0.000076, 0.017, 0.07, 0.145, 0.245, 0.371, 0.5, 0.627, 0.75, 0.852, 0.932, 0.983, 1]

	total_count = desired_number / interpolate(angle, xs, ys)
	angle2 = abs(angle)

	vector_list = make_vectors_sphere(int(total_count))
	selected_vectors = [v for v in vector_list if abs(calculate_angle(vector, v, True)) < angle2]

	return selected_vectors


def get_PCA_vector(helices):
	vectors = []
	for helix_vector in helices:
		vectors.append(np.array(helix_vector))
		vectors.append(-np.array(helix_vector))
	pca1 = PCA()
	pca1.fit(vectors)
	return pca1.components_[0]
