from pymol import cmd
from pymol.vfont import plain
from pymol.cgo import *

import numpy as np
import json


def make_line(p1, p2, color=[1,1,1], thickness=0.3):
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    return [CYLINDER, float(x1), float(y1), float(z1), 
    float(x2), float(y2), float(z2), 
    thickness, color[0], color[1], color[2], color[0], color[1], color[2]]

def make_normal_line(normal, position, thickness, color1=[1,1,1]):
	npoint1 = position + thickness * 1.5 * normal
	npoint2 = position - thickness * 1.5 * normal
	lineObj = []
	lineObj.extend(make_line(npoint1, npoint2, color1, 0.5))
	return lineObj


def get_plane_vectors(normal):
	"""
	From a normal vector, computes 2 orthogonal vectors v1 and v2 (normed)
	for spanning the plane of the normal vector.
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

def make_plane(corner1, corner2, corner3, corner4, normal, color=[0.8, 0.8, 0.8], alpha=0.6):
	planeObj = []
	planeObj.extend([ALPHA, alpha])
	planeObj.extend([COLOR, color[0], color[1], color[2]])

	planeObj.extend([BEGIN, TRIANGLE_STRIP])
	planeObj.append(NORMAL)
	planeObj.extend(normal)

	for corner in [corner1, corner2, corner3, corner4, corner1]:
		planeObj.append(VERTEX)
		planeObj.extend(corner)
	planeObj.append(END)

	return planeObj

def make_membrane(normal, position, thickness, size, name, color=[0.8, 0.8, 0.8], alpha=0.6):
	v1, v2 = get_plane_vectors(normal)

	# Make Corner points of planes
	corner1_bot = position - thickness * normal + v1 * size
	corner2_bot = position - thickness * normal + v2 * size
	corner3_bot = position - thickness * normal - v1 * size
	corner4_bot = position - thickness * normal - v2 * size

	corner1_top = position + thickness * normal + v1 * size
	corner2_top = position + thickness * normal + v2 * size
	corner3_top = position + thickness * normal - v1 * size
	corner4_top = position + thickness * normal - v2 * size

	# Make Plane Object
	plane1 = make_plane(corner1_bot, corner2_bot, corner3_bot, corner4_bot, normal, color, alpha)
	plane2 = make_plane(corner1_top, corner2_top, corner3_top, corner4_top, normal, color, alpha)
	normal_line = make_normal_line(normal, position, thickness, color)

	makePymolObject(normal_line, "normal_"+name)
	makePymolObject(plane1, "plane1_"+name)
	makePymolObject(plane2, "plane2_"+name)
	cmd.group(name, "normal_"+name+" plane1_"+name+" plane2_"+name)


def makePymolObject(cgo, name):
	auto_zoom = cmd.get('auto_zoom', quiet=1)
	cmd.set('auto_zoom', 0, quiet=1)
	cmd.load_cgo(cgo, name)
	cmd.set('auto_zoom', auto_zoom, quiet=1)




def visualise(path_to_pdb_file, path_to_json):
	# Load everything
	cmd.reinitialize()
	cmd.load(path_to_pdb_file)
	pdb_id = path_to_pdb_file.split('/')[-1].split('.')[0]

	with open(path_to_json) as file:
	    data = json.load(file)

	# Position/visualise PDB model
	cmd.translate([-coord for coord in data['init_pos']], pdb_id)
	cmd.orient(pdb_id)
	cmd.show("cartoon", pdb_id)
	cmd.hide("lines", pdb_id)
	cmd.color("red", "ss h")
	cmd.color("yellow", "ss s")
	cmd.set("cartoon_transparency", 0.6, pdb_id)


	# All helices and TM helices
	obj_all = []
	obj_tm = []

	for v, c, l in zip(data['vectors_all'], data['centers_all'], data['labels_all']):
		v2 = np.array(v)/np.linalg.norm(np.array(v))*15
		c = np.array(c)

		obj_all.extend(make_line(c - v2, c + v2))
		if l==1:
			obj_tm.extend(make_line(c - v2, c + v2, [0.3,0.3,1]))

	makePymolObject(obj_all, "helices_all")
	makePymolObject(obj_tm, "helices_TM")


	# Filter Step 1: by angle
	pca1_p1 = np.array(data['pca1'])*25
	pca1_p2 = -pca1_p1
	obj_pca1 = make_line(pca1_p1, pca1_p2, color=[0.5, 0.5, 1], thickness=1)
	makePymolObject(obj_pca1, "normal_pre_filter")
	
	obj_filter1 = []
	for v, c in zip(data['vectors_filter1'], data['centers_filter1']):
		v2 = np.array(v)/np.linalg.norm(np.array(v))*15
		c = np.array(c)
		obj_filter1.extend(make_line(c - v2, c + v2, color=[0,1,0]))

	makePymolObject(obj_filter1, "helices_filter1")


	# Filter Step 2: by projected distance
	obj_filter2 = []
	for v, c in zip(data['vectors_filter2'], data['centers_filter2']):
		v2 = np.array(v)/np.linalg.norm(np.array(v))*15
		c = np.array(c)
		obj_filter2.extend(make_line(c - v2, c + v2, color=[1,1,0]))

	makePymolObject(obj_filter2, "helices_filter2")


	# Membranes
	thickness = 17
	size = 100
	# Membrane 1: Before Optimisation - average of TM helices
	normal1 = np.array(data['vector_final'])
	position1 = np.array(data['center_final'])
	
	make_membrane(normal1, position1, thickness, size, 'prelim', color=[1, 0.3, 1], alpha=0.6)

	# Membrane 2: After Optimisation
	normal2 = np.array(data['vector_opti'])
	position2 = data['z_opti']*normal2
	make_membrane(normal2, position2, thickness, size, 'opti', color=[0.3, 1, 0.3], alpha=0.6)


	cmd.disable("helices_all helices_filter1 normal_pre_filter")


print("Functions loaded. Use   visualise('path_to_pdb_file', 'path_to_json')   to visualise a membrane.")

#cmd.show("cgo", "plane*")
