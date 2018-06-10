from pymol import cmd
from pymol.vfont import plain
from pymol.cgo import *
import numpy as np

cmd.delete("atom1 AND atom2 AND atom3 AND atom4 AND plane* AND normal-*")

def line(p1, p2, color1=[1,1,1], color2=[1,1,1]):
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    return [CYLINDER, float(x1), float(y1), float(z1), float(x2), float(y2), float(z2), 0.25, 
    color1[0], color1[1], color1[2], color2[0], color2[1], color2[2]]

def get_plane_vectors(normal):
	"""
	From a normal vector, computes 2 orthogonal vectors
	for spanning the plane of the normal vector
	:param normal: 3D numpy array
	:return: v1 - 3D numpy array
	         v2 - 3D numpy array
	"""
	if normal[1] == 0 and normal[2] == 0:
	    if normal[0] == 0:
	        raise ValueError('Normal vector is (0,0,0) -> no plane possible')
	    else:
	        v1 = np.cross(normal, [0, 1, 0])
	        v2 = np.cross(normal, v1)
	        return v1, v2
	v1 = np.cross(normal, [1, 0, 0])
	v2 = np.cross(normal, v1)
	return v1/np.linalg.norm(v1), v2/np.linalg.norm(v2)


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


	planeObj.extend(-1*normal)

	for corner in [corner1, corner2, corner3, corner4, corner1]:
		planeObj.append(VERTEX)
		planeObj.extend(corner)
	planeObj.append(END)

	return planeObj


def makeLine(normal, position, thickness, color1=[1,1,1], color2=[1,1,1]):
	npoint1 = position + thickness * 1.5 * normal
	npoint2 = position - thickness * 1.5 * normal
	lineObj = []
	lineObj.extend(line(npoint1, npoint2, color1, color2))
	
	return lineObj



def makePrimitive(cgo, name):
	az = cmd.get('auto_zoom', quiet=1)
	cmd.set('auto_zoom', 0, quiet=1)
	cmd.load_cgo(cgo, name)
	cmd.set('auto_zoom', az, quiet=1)


def main():
	thickness = 17
	size = 80
	cmd.reinitialize()
	cmd.load("1be3.pdb")
	cmd.translate([ -32.23923999, -127.88527506,  -21.56395075], "1be3")
	cmd.pseudoatom("ori", pos=[0,0,0])
	cmd.show("cartoon", "1be3")
	cmd.hide("lines", "1be3")
	cmd.orient("1be3")
	cmd.color("red", "ss h")
	cmd.color("yellow", "ss s")

	# Averaging the TM helices

	position = np.array( [ 15.18373757, -8.18554451, 13.23525919])
	normal = np.array([0.79483971, -0.60645899,  0.02091262])

	v1, v2 = get_plane_vectors(normal)

	corner1_bot = position - thickness * normal + v1 * size
	corner2_bot = position - thickness * normal + v2 * size
	corner3_bot = position - thickness * normal - v1 * size
	corner4_bot = position - thickness * normal - v2 * size

	corner1_top = position + thickness * normal + v1 * size
	corner2_top = position + thickness * normal + v2 * size
	corner3_top = position + thickness * normal - v1 * size
	corner4_top = position + thickness * normal - v2 * size
	settings = {'ALPHA': 0.6, 'COLOR': [1, 0.3, 1]}

	plane1 = plane(corner1_bot, corner2_bot, corner3_bot, corner4_bot, normal, settings)
	plane2 = plane(corner1_top, corner2_top, corner3_top, corner4_top, normal, settings)
	normal_line = makeLine(normal, position, thickness, [1, 0.3, 1], [1, 0.3, 1])
	makePrimitive(normal_line, "normal-ave")
	makePrimitive(plane1, "plane1-ave")
	makePrimitive(plane2, "plane2-ave")


		
	#position = np.array([45.920853, 106.799995,  35.801376])
	normal = np.array([0.7034874344833971, -0.6717687074829933, 0.23201774322365207])
	position = 16*normal

	v1, v2 = get_plane_vectors(normal)

	corner1_bot = position - thickness * normal + v1 * size
	corner2_bot = position - thickness * normal + v2 * size
	corner3_bot = position - thickness * normal - v1 * size
	corner4_bot = position - thickness * normal - v2 * size

	corner1_top = position + thickness * normal + v1 * size
	corner2_top = position + thickness * normal + v2 * size
	corner3_top = position + thickness * normal - v1 * size
	corner4_top = position + thickness * normal - v2 * size
	settings = {'ALPHA': 0.6, 'COLOR': [0.3, 1, 0.3]}

	plane1 = plane(corner1_bot, corner2_bot, corner3_bot, corner4_bot, normal, settings)
	plane2 = plane(corner1_top, corner2_top, corner3_top, corner4_top, normal, settings)
	normal_line = makeLine(normal, position, thickness, [0.3, 1, 0.3], [0.3, 1, 0.3])
	makePrimitive(normal_line, "normal-opt")
	makePrimitive(plane1, "plane1-opt")
	makePrimitive(plane2, "plane2-opt")



	cmd.show("cgo", "normal-*")
	cmd.show("cgo", "plane*")
	cmd.group("opti", "normal-opt plane1-opt plane2-opt")
	cmd.group("ave", "normal-ave plane1-ave plane2-ave")



main()
