#!/home/mat/.pyenv/shims/python python3
# system imports
import sys
import os
import math
import time

# required packages
import regex as re
import numpy as np
import mcubes

# local imports
# these require numpy
import scalar_field as sf
import molecule as mol

class CubeSettings():
	def __init__(self):
		self.roll = False

class Cube():
	""" CLASS Cube()
	cube object stores everything in the cube file, to be loaded in as efficiently as possible
	that is, for large objects, loading is done lazily
	"""

	def __init__(self):
		# all information relating to the field is stored here
		self.field = sf.ScalarField()
		# all atomic information is stored here
		self.molecule = mol.Molecule()
		self.settings = CubeSettings()
		self.file = None
		self.name = None

	def __str__(self):
		ostr = "FIELD DETAILS:\n"
		ostr += self.field.__str__() + "\n"
		ostr += "MOLECULAR DETAILS:\n"
		ostr += self.molecule.__str__()
		return ostr

	def load_header(self, file):
		""" FUNCTION load_header(string: file)
		reads a cube input (if not cube, it will complain) and only extracts header, along with pointers to
		the file location.
		"""
		abspath = os.path.abspath(file)
		self.file = abspath
		self.name = abspath.split('/')[-1].split('.')[0]
		self.field.load_file(abspath)
		self.molecule.load_file(abspath)

		with open(file, 'r') as f:
			atomcount = 0
			for idx, line in enumerate(f):
				if idx < 2:
					# This is currently the expected output from Quantum ESPRESSO with Environ
					# There may be more or less text here, so TODO generalize
					continue
				if idx == 2:
					# this line contains the number of atoms and the position of the origin
					# split by 1 or more spaces
					l = re.split('\s+', line.strip())
					self.molecule.load_empty(l[0])
					atomcount = int(l[0])
					self.field.set_translation(l[1:])
					continue
				if idx < 6:
					# the next 3 lines contain the size and transformation matrices
					l = re.split('\s+', line.strip())
					self.field.add_sizeparam(idx-3, l[0])
					self.field.add_transform(idx-3, l[1:])
					continue
				if idx < (6 + atomcount):
					# the next 'atomcount' lines contain the atomic information
					l = re.split('\s+', line.strip())
					self.molecule.add_atom(idx-6, l)
					continue
				if idx > (5 + atomcount):
					break

		self.field.add_scaling()
		self.molecule.transform(self.field.transform)

	def field_settings(self, roll=False):
		self.settings.roll = roll

	def load_body(self):
		"""
		reads in the field
		"""
		start = time.time()
		print('reading field...')
		self.field.init_field()
		with open(self.file, 'r') as f:
			fx = 0
			fy = 0
			fz = 0
			for idx, line in enumerate(f):
				if idx < (6 + self.molecule.atomcount):
					continue
				l = re.split('\s+', line.strip())
				for e in l:
					if e is None or e == '':
						continue
					elif e == 'NaN':
						self.field.field[fx, fy, fz] = 0.0
					else:
						self.field.field[fx, fy, fz] = e
					fx, fy, fz = self.field.increment_idx(fx, fy, fz)
					if fx < 0:
						break

		if self.settings.roll:
			self.field.roll()

		end = time.time()
		print('field loading complete, time elapsed = {}s'.format(end-start))

	def save_voxel(self, path, voxeldata):
		# create header
		header = np.zeros((4,), dtype=int)
		header[0:3] = self.field.gridsize
		header[3] = 1 # for still frame

		with open(path, 'wb') as binfile:
			header.astype('<i4').tofile(binfile)
			voxeldata.astype('<f4').tofile(binfile)

	def check_file(self, name, update):
		current_dir = os.path.dirname(os.path.realpath(__file__))
		# if dat folder doesn't exist, something is wrong but create it anyway
		if not os.path.isdir(os.path.join(current_dir, 'dat')):
			os.mkdir('dat')
		# check dat folder for occurances
		vpath = os.path.join(current_dir, 'dat', name)
		if not update and os.path.isfile(vpath):
			print('{} already exists'.format(name))
			return None # nothing needed to be done
		print('{} does not exist'.format(vpath))
		if self.field.field is None:
			self.load_body()
		return vpath

	def make_isomesh(self, val, name="", update=False):
		if not name:
			name = self.name + '.dae'
		elif not name.endswith('.dae'):
			name = name + '.dae'
		ipath = self.check_file(name, update)
		if ipath is None:
			return
		print('making isosurface...')
		start = time.time()
		field_max = np.amax(self.field.field)
		field_min = np.amin(self.field.field)
		isoval = val * (field_max - field_min) + field_min
		vertices, triangles = mcubes.marching_cubes(self.field.field, isoval)
		mcubes.export_mesh(vertices, triangles, ipath, "Iso{}".format(val))
		end = time.time()
		print('mesh created, time elapsed = {}s'.format(end-start))

	# creating isosurfaces and voxel files are expensive. Save the files for repeat use.
	def make_color_voxel(self, name="", update=False):
		if not name:
			name = self.name + '_color.bvox'
		elif not name.endswith('.bvox'):
			name = name + '_color.bvox'
		vpath = self.check_file(name, update)
		if vpath is None:
			return
		print('making color voxel...')
		start = time.time()
		vox = self.field.field.flatten()
		# normalize
		vox -= np.min(vox)
		vox /= np.max(vox)
		# flip
		#vox = 1.0 - vox
		# save
		self.save_voxel(vpath, vox)
		end = time.time()
		print('color voxel created, time elapsed = {}s'.format(end-start))

	def make_emission_voxel(self, name="", update=False, truncA=-1e20, truncB=1e20,
							max_emission=0.5, tol=0.1, modifier='SIGMOID'):
		if not name:
			name = self.name + '_emission.bvox'
		elif not name.endswith('.bvox'):
			name = name + '_emission.bvox'
		vpath = self.check_file(name, update)
		if vpath is None:
			return
		print('making emission voxel...')
		start = time.time()
		field = self.field.field
		# clip if desirable
		np.clip(field, truncA, truncB)
		# normalize
		field -= np.min(field)
		field /= np.max(field)
		# flip
		field = 1.0 - field
		if modifier == 'SIGMOID':
			# convert to sigmoid input
			field -= 0.5
			field *= 8.0
			# apply sigmoid function
			field = np.exp(field)
			field /= (field + 1)
			field = np.clip(field, 0, max_emission)
			field[field < tol] = 0.0
		elif modifier == 'GRADIENT':
			print('non-standard gradient modifier chosen')
			# to save computation time, just do a rough finite forward difference on the three adjacent
			# cells (x, y, z)
			finite_x = (np.roll(field, 1, axis=0) + np.roll(field, 1, axis=1) + np.roll(field, 1, axis=2)
						- (3 * field)) / 3.0
			field = np.clip(field, 0, max_emission)
		else:
			# not recognized, just do nothing and hope for the best...
			print('warning: modifier option not recognized')
		# save
		vox = field.flatten()
		self.save_voxel(vpath, vox)
		end = time.time()
		print('emission voxel created, time elapsed = {}s'.format(end-start))

if __name__ == '__main__':
	cube = Cube()
	cube.load_header('test.cube')
	print(cube)
	cube.load_body()
	cube.make_color_voxel(update=True)
	cube.make_emission_voxel(update=True)
