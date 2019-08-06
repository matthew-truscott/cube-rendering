#!/home/mat/.pyenv/shims/python python3

import numpy as np

class ScalarField():
	"""
	scalar field object that contains transformation details from grid space to real space, and
	a single floating point value for each gridspce.

	Initialization is skeletal (that is, the field contents do not get loaded in and instead, a
	pointer to the cube file is stored until the field needs to be read in).
	"""

	def __init__(self):
		self.status = 0 # how much of the field is initialized
		self.gridsize = np.zeros((3,), dtype=int)
		self.field = None
		self.transform = np.eye(4, dtype=float)
		self.meshtransform = None
		self.cubefile = None

	def __str__(self):
		ostr = '\tGRIDSIZE = {}'.format(self.gridsize) + '\n'
		ostr += '\tFIELDPATH = {}'.format(self.cubefile) + '\n'
		ostr += '\tTRANSFORM = \n{}'.format(self.transform) + '\n'
		return ostr

	def roll(self):
		# roll the field
		uroll = self.gridsize // 2
		self.field = np.roll(self.field, uroll, axis=(0, 1, 2))

	def init_field(self):
		self.field = np.zeros((self.gridsize), dtype=float)

	def load_file(self, file):
		self.cubefile = file

	def set_translation(self, translate):
		self.transform[3, 0:3] = translate

	def add_sizeparam(self, idx, val):
		self.gridsize[idx] = val

	def add_transform(self, idx, vallist):
		self.transform[idx, 0:3] = vallist

	def add_scaling(self):
		scale = np.zeros((4,), dtype=int)
		scale[0:3] = self.gridsize
		self.meshtransform = np.copy(self.transform)
		self.transform = self.transform * scale[:, np.newaxis]

	def increment_idx(self, x, y, z):
		z += 1
		if z >= self.gridsize[2]:
			z = 0
			y += 1
			if y >= self.gridsize[1]:
				y = 0
				x += 1
				if x >= self.gridsize[0]:
					return -1, -1, -1
		return x, y, z
