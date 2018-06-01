import math


def make_vectors(n=100):
	"""
	Generates an array of evenly distributed unit vectors
	in 3D Space	using a Fibonacci sphere.
	Adapted from
		https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere/26127012#26127012
	:param n: Number of vectors to be generated
	:return: List of generated vectors
	"""
	vectors = []
	y_offset = 2.0 / n
	step_size = math.pi * (3.0 - math.sqrt(5.0))

	for i in range(n):
		y = ((i * y_offset) - 1) + (y_offset / 2)
		r = math.sqrt(1 - y*y)

		phi = ((i + 1.0) % n) * step_size

		x = math.cos(phi) * r
		z = math.sin(phi) * r

		vectors.append([x,y,z])

	return vectors
