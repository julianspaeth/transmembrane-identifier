from pymol import cmd
from pymol.vfont import plain
from pymol.cgo import *
import numpy as np

cmd.delete("atom1 AND atom2 AND atom3 AND atom4")


settings = {'ALPHA': 0.6, 'COLOR': [1, 1, 0], 'INVERT': True}

def plane(corner1, corner2, corner3, corner4, normal, settings):
	planeObj = []

	if 'ALPHA' in settings:
		planeObj.extend([ALPHA, settings['ALPHA']])

	if 'COLOR' in settings:
		planeObj.extend([COLOR, settings['COLOR'][0], settings['COLOR'][1], settings['COLOR'][2]])
	else:
		planeObj.extend([COLOR, 0.8, 0.8, 0.8]) # greyish

	planeObj.extend([BEGIN, TRIANGLE_STRIP])
	planeObj.append(NORMAL)

	if 'INVERT' in settings:
		if settings['INVERT']==True:
			planeObj.extend(-1*normal)
		else:
			planeObj.extend(normal)
	else:
		planeObj.extend(normal)

	for corner in [corner1, corner2, corner3, corner4, corner1]:
		planeObj.append(VERTEX)
		planeObj.extend(corner)
	planeObj.append(END)
	return planeObj


def planeFromPoints(point1, point2, point3, size):
	v1 = np.subtract(point2, point1)
	v1 = v1 / np.linalg.norm(v1)

	v2 = np.subtract(point3, point1)
	v2 = v2 / np.linalg.norm(v2)

	normal = np.cross(v1, v2)

	v2_new = np.cross(normal, v1)

	x = v1 * size
	y = v2 * size
	center = point2

	corner1 = np.add(np.add(center, x), y)
	corner2 = np.subtract(np.add(center, x), y)
	corner3 = np.subtract(np.subtract(center, x), y)
	corner4 = np.add(np.subtract(center, x), y)
	return plane(corner1, corner2, corner3, corner4, normal, dicti)


def makePrimitive(cgo, name):
	az = cmd.get('auto_zoom', quiet=1)
	cmd.set('auto_zoom', 0, quiet=1)
	cmd.load_cgo(cgo, name)
	cmd.set('auto_zoom', az, quiet=1)


def main():
	size = 10
	plane = planeFromPoints(point1, point2, point3, size)
	planeName = "plane-1"
	makePrimitive(plane, planeName)
	cmd.show("cgo", "plane*")


point1 = np.array([0,0,0])
point2 = np.array([10,0,0])
point3 = np.array([10,10,0])

main()
