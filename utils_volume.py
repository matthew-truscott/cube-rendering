import sys
import os
import bpy
import cube_reader as cr
from math import pi

def add_isosurface(cube, val, name="", update=False):
	""" FUNCTION add_isosurface(cube: Cube, name: str, update: bool)
	Adds an isosurface object using the marching cubes external package (may want to implement this
	in the future for more freedom, but it works fine for now).

	INPUT:
	Cube: cube, the object that contains all the cell data
	float: val, a value between 0 and 1 that will determine the field. TODO absolute val option here too?

	RETURNS:
	blender object: the isosurface that represents the field data from the relevant cube file, taken
	at a particular value
	"""
	if not name:
		name = cube.name + '.dae'
	else:
		name = name + '.dae'
	# assume mesh file does not exist, so run external checker. If they indeed do exist, try
	# loading the relevant cube file referenced by 'name' and create the files on the fly
	cube.make_isomesh(val, name, update=update)

	# set directories
	current_dir = os.path.dirname(os.path.realpath(__file__))
	data_dir = os.path.join(current_dir, 'dat')
	isodir = os.path.join(data_dir, name)

	# cube position needs to be fixed for the mesh
	cube_position = cube.field.gridsize / 2.0
	cube_position[1] = -cube_position[1]
	print('cubeposition', cube_position)
	bpy.context.scene.cursor_location = cube_position
	bpy.ops.wm.collada_import(filepath=isodir)
	obj = bpy.context.active_object
	bpy.ops.object.origin_set(type='ORIGIN_CURSOR')

	# apply transforms to scale appropriately to the molecule drawing
	obj.rotation_euler = (0, 0, 0)
	obj.data.transform(cube.field.meshtransform)
	obj.location = [0, 0, 0]

	# default look for the surface..
	mat = bpy.data.materials.new('SurfaceMaterial')
	mat.use_transparency = True
	mat.transparency_method = 'RAYTRACE'
	mat.raytrace_transparency.fresnel = 4.0
	obj.data.materials.append(mat)
	obj.data.update()

	return obj

def add_volume(cube, name="", update=False):
	""" FUNCTION add_volume(cube: Cube, name: str)
	Adds a volume object for blender to render, based off the voxel data from a cube file
	Requires the existence of these voxel files, which are handled by another function (see
	cube_reader.py)

	INPUT:
	Cube: cube, the object that contains all the cell data
	str: name, the name of the cube/voxel file to be read (these should have the same name, but if not,
		reference the voxel file name(s, as a list)

	RETURNS:
	blender object: the cube that represents the field data from the relevant cube file, referenced by
	the name given.
	"""
	if not name:
		# default expectation
		name0 = cube.name + '_color.bvox'
		name1 = cube.name + '_emission.bvox'
		cube.make_color_voxel(update=update)
		cube.make_emission_voxel(update=update, modifier='GRADIENT', max_emission=0.05)
	elif isinstance(name, (list,)):
		if len(name) != 2:
			print("name parameter expects 2-list or string")
			return
		else:
			if name[0].endswith('.bvox'):
				name0 = name[0]
			else:
				# try adding extension
				name0 = name[0] + '.bvox'
			if name[1].endswith('bvox'):
				name1 = name[1]
			else:
				name1 = name[1] + '.bvox'
	else:
		# assume names are appended with color and emission
		name0 = name + '_color.bvox'
		name1 = name + '_emission.bvox'
		# assume voxel files do not exist, so run external checker. If they indeed do not exist, 
		# try loading the relevant cube file referenced by 'name' and create the files on the fly.
		cube.make_color_voxel(name, update=update)
		cube.make_emission_voxel(name, update=update, modifier='GRADIENT', max_emission=0.2)

	# set directories
	current_dir = os.path.dirname(os.path.realpath(__file__))
	data_dir = os.path.join(current_dir, 'dat')

	# add primitive cube and set as main object for convenient manipulation
	bpy.ops.mesh.primitive_cube_add(location=(0, 0, 0))
	obj = bpy.context.active_object

	# set references to useful objects in cube
	cellsize = cube.field.gridsize
	cellt = cube.field.transform

	# fortran stores field in inverse order such that x and z are flipped when read directly by
	# blender, the following two transformations fixes this, along with a rescaling since the blender
	# defaults to length 2.0
	obj.scale = (-0.5, 0.5, 0.5)
	obj.rotation_euler = (0, pi / 2.0, 0)
	obj.data.transform(cellt)
	obj.data.update()

	# create material based off volume
	mat = bpy.data.materials.new('VolumeMaterial')
	mat.type = 'VOLUME'
	mat.volume.scattering = 0.05
	mat.volume.reflection = 0.0

	# texture 0 sets the color of the material, based off values of voxel data
	tex0 = bpy.data.textures.new('VoxelTexture', 'VOXEL_DATA')
	tex0.voxel_data.file_format = 'BLENDER_VOXEL'
	tex0.voxel_data.filepath = os.path.join(data_dir, name0)
	tex0.use_color_ramp = True

	# texture 1 sets the emission of the material, based off a scaled dataset, for a more defined visible
	# volume
	tex1 = bpy.data.textures.new('VoxelTexture', 'VOXEL_DATA')
	tex1.voxel_data.file_format = 'BLENDER_VOXEL'
	tex1.voxel_data.filepath = os.path.join(data_dir, name1)

	# add texture 0 to the material
	slot0 = mat.texture_slots.add()
	slot0.texture = tex0
	slot0.texture_coords = 'ORCO'
	slot0.use_map_color_emission = True

	# add texture 1 to the material
	slot1 = mat.texture_slots.add()
	slot1.texture = tex1
	slot1.texture_coords = 'ORCO'
	slot1.use_map_emission = True

	obj.data.materials.append(mat)

	# return object, everything else should be able to be referenced from this object
	return obj

def set_volume_color(obj):
	mat = obj.data.materials[0]
	tex = mat.texture_slots[0].texture # 0th slot should be reserved for color

	# color map, should be defined externally
	cr0 = tex.color_ramp.elements.new(0.0)
	cr1 = tex.color_ramp.elements.new(0.3)
	cr2 = tex.color_ramp.elements.new(1.0)

	cr0.color = (0, 0, 0, 0)
	cr1.color = (0, 0, 10, 0.1)
	cr2.color = (0, 0, 10, 0.1)
