import bpy
import sys
import os
import numpy as np

dir = os.path.dirname(bpy.data.filepath)
if not dir in sys.path:
	sys.path.append(dir)

print(sys.path)
print(sys.executable)

import utils_blender as ub
import utils_volume as uv
import utils_molecule as um
import cube_reader as cr

def create_scene():
	# add the molecule
	cube = cr.Cube()
	# load the molecule
	cube.load_header('cube/3.cube')
	#cube.field_settings(roll=True)

	# clean up
	ub.removeAll()

	# set the scene
	target = ub.target()
	scene = ub.get_scene()
	camera = ub.camera((0, 0, 30), target)
	scene.camera = camera
	lamp = ub.lamp((-14, -16, 18), type='POINT', energy=4, color=(1,1,1), target=target)
	lamp = ub.lamp((13, 15, -10), type='POINT', energy=4, color=(1,1,1), target=target)
	bpy.data.worlds['World'].color = (0, 0, 0)

	# TODO make volume a class so that one can successively render objects with persistence
	#vol = uv.add_volume(cube, update=False)
	#uv.set_volume_color(vol)
	#iso = uv.add_isosurface(cube, 0.5, update=True)
	#iso = uv.add_isosurface(cube, 0.6, update=True)
	# for the above (more than one call) need to have smarter filenaming, and
	# need to change the max allowed reflections for reasonable transparency

	# update molecule in order to get pointers to names of rendered objects, for future editing
	cube = um.draw_molecule(cube, bonds=True)
	um.edit_atom_material(cube, specular=0.4)

	ub.save('blendertest')
	ub.renderToFolder('test', '01', 800, 800)

if __name__ == '__main__':
	create_scene()
