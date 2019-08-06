#!/home/mat/.pyenv/shims/python python3

import numpy as np

class Molecule():
	"""
	molecule contains all the molecular information along with functions that can be called to aid with
	molecular rendering
	"""

	def __init__(self):
		self.status = 0 # how much of the molecule has been initialized
		self.atomcount = 0
		self.a_species = None
		self.bondlist = []
		self.a_charges = None
		self.m_positions = None
		self.label = "" # optional label to be set
		self.cubefile = None
		self.rendered = []

	def __str__(self):
		ostr = '\tATOMCOUNT = {}\n'.format(self.atomcount)
		ostr += '\tSPECIES LIST = {}\n'.format(self.a_species)
		ostr += '\tCHARGES LIST = {}\n'.format(self.a_charges)
		ostr += '\tPOSITIONS =\n{}'.format(self.m_positions)
		return ostr

	def load_empty(self, atomcount):
		# return status error if atomcount is invalid
		try:
			atomcount = int(atomcount)
		except ValueError:
			return 1
		self.status = 1
		self.atomcount = atomcount
		self.a_species = np.zeros((atomcount), dtype=int)
		self.a_charges = np.zeros((atomcount), dtype=float)
		self.m_positions = np.zeros((atomcount, 3), dtype=float)

		return 0

	def load_file(self, file):
		self.cubefile = file

	def add_atom(self, idx, datalist):
		self.a_species[idx] = datalist[0]
		self.a_charges[idx] = datalist[1]
		self.m_positions[idx] = datalist[2:5]

	def center_molecule(self):
		center = np.mean(self.m_positions, axis=0)
		self.m_positions = self.m_positions - center

	def transform(self, transform):
		self.m_positions = self.m_positions - (transform.diagonal()/2)[0:3]

	def get_distance(self, vecx, vecy):
		sdx = (vecx[0] - vecy[0]) ** 2
		sdy = (vecx[1] - vecy[1]) ** 2
		sdz = (vecx[2] - vecy[2]) ** 2
		dist = (sdx + sdy + sdz) ** 0.5
		pmdist = dist * 52.91772083
		return pmdist

	def create_bonds(self):
		# have lists for the atomic positions and the atomlist. 
		# Now, connect the atoms using expected bond lengths
		# Each atom has a valence, may want to figure this out 
		# in the future (but this isn't a general rule because
		# ions
		# Instead just loop through entire list each time

		# bond lists for common bonding atoms, index is equivalent to atomic mass, each entry should be a tuple of
		# min/max in pm
		bond_dict = {}
		bond_dict[6] = {}
		bond_dict[6][1] = (106, 112)
		bond_dict[6][6] = (120, 154)
		bond_dict[6][7] = (116, 210)
		bond_dict[6][8] = (113, 215)
		bond_dict[7] = {}
		bond_dict[7][1] = (90, 110)
		bond_dict[7][7] = (120, 155)
		bond_dict[8] = {}
		bond_dict[8][1] = (90, 100)
		bond_dict[8][6] = (120, 154)

		self.bondlist = []
		atomcheck = set(range(self.atomcount))
		# iterate through the list of positions, get distance between atoms and compare with dicts
		for i, ival in enumerate(self.a_species):
			for j, jval in enumerate(self.a_species):
				# skip repeats
				if j <= i:
					continue
				min_length = 0
				max_length = 0
				if ival in bond_dict:
					if jval in bond_dict[ival]:
						min_length, max_length = bond_dict[ival][jval]
				elif jval in bond_dict:
					if ival in bond_dict[jval]:
						min_length, max_length = bond_dict[jval][ival]
				distance = self.get_distance(self.m_positions[i], self.m_positions[j])
				#print('distance between {} and {} is {}'.format(i, j, distance))
				if distance > min_length and distance < max_length:
					self.bondlist.append((i, j))
					#print('bonded {} with {}'.format(i, j))
					if i in atomcheck:
						atomcheck.remove(i)
					if j in atomcheck:
						atomcheck.remove(j)

		print('bonding complete: {}'.format(self.bondlist))
		if atomcheck:
			print('the following indices are unattached: {}'.format(atomcheck))

	def add_rendered(self, atoms):
		self.rendered = atoms
