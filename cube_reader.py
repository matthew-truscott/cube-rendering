"""CUBE_READER

Author: Matthew Truscott
"""
#!/home/mat/.pyenv/shims/python python3

# system imports
import os
import time
import sys

# required packages
import mcubes
import regex as re
import numpy as np

# local imports
import scalar_field as sf
import molecule as mol

class CubeSettings():
    """cube reading settings container
    """
    def __init__(self):
        self.roll = False

class Cube():
    """.cube object container
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
        reads a cube input (if not cube, it will complain) and only extracts header, along with 
        pointers to the file location.
        """
        abspath = os.path.abspath(file)
        self.file = abspath
        self.name = abspath.split('/')[-1].split('.')[0]
        self.field.load_file(abspath)
        self.molecule.load_file(abspath)

        with open(file, 'r') as f_read:
            atomcount = 0
            for idx, line in enumerate(f_read):
                if idx < 2:
                    # This is currently the expected output from Quantum ESPRESSO with Environ
                    # There may be more or less text here, so TODO generalize
                    continue
                if idx == 2:
                    # this line contains the number of atoms and the position of the origin
                    # split by 1 or more spaces
                    line_elements = re.split(r'\s+', line.strip())
                    self.molecule.load_empty(line_elements[0])
                    atomcount = int(line_elements[0])
                    self.field.set_translation(line_elements[1:])
                    continue
                if idx < 6:
                    # the next 3 lines contain the size and transformation matrices
                    line_elements = re.split(r'\s+', line.strip())
                    self.field.add_sizeparam(idx-3, line_elements[0])
                    self.field.add_transform(idx-3, line_elements[1:])
                    continue
                if idx < (6 + atomcount):
                    # the next 'atomcount' lines contain the atomic information
                    line_elements = re.split(r'\s+', line.strip())
                    self.molecule.add_atom(idx-6, line_elements)
                    continue
                if idx > (5 + atomcount):
                    break

        self.field.add_scaling()
        self.molecule.transform(self.field.transform)


    def field_settings(self, roll=False):
        """Settings for the field container
        
        Keyword Arguments:
            roll {bool} -- Roll the field in order to correctly display isolated
            molecules that are defined across cell edges (default: {False})
        """
        self.settings.roll = roll


    def load_body(self):
        """Reads in the field
        """
        start = time.time()
        print('reading field...')
        self.field.init_field()
        with open(self.file, 'r') as f_read:
            f_x = 0
            f_y = 0
            f_z = 0
            for idx, line in enumerate(f_read):
                if idx < (6 + self.molecule.atomcount):
                    continue
                line_elements = re.split(r'\s+', line.strip())
                for element in line_elements:
                    if element is None or element == '':
                        continue
                    elif element == 'NaN':
                        self.field.field[f_x, f_y, f_z] = 0.0
                    else:
                        self.field.field[f_x, f_y, f_z] = element
                    f_x, f_y, f_z = self.field.increment_idx(f_x, f_y, f_z)
                    if f_x < 0:
                        break

        if self.settings.roll:
            self.field.roll()

        end = time.time()
        print('field loading complete, time elapsed = {}s'.format(end-start))


    def save_voxel(self, path, voxeldata):
        """saves the data to a voxel file
        
        Arguments:
            path {string} -- path of the save file
            voxeldata {np array} -- the data for the voxel file
        """
        # create header
        header = np.zeros((4,), dtype=int)
        header[0:3] = self.field.gridsize
        header[3] = 1 # for still frame

        with open(path, 'wb') as binfile:
            header.astype('<i4').tofile(binfile)
            voxeldata.astype('<f4').tofile(binfile)


    def check_file(self, name, update):
        """Checks that the file exists, making this data file is expensive and
        the file is large, so avoid if possible!

        Arguments:
            name {string} -- name of the file to be checked
            update {bool} -- update if exists?

        Returns:
            path -- the path of the file
        """
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
        """makes a mesh based off the marching cubes algorithm, for given volume data
        
        Arguments:
            val {float} -- value between 0 and 1 that determines where the mesh is drawn. Mesh is
            an isosurface based off volumetric data. 0 takes the minimum value in the volume and
            tries to make a surface on that value, 1 takes the maximum. 
        
        Keyword Arguments:
            name {str} -- given name for isomesh (default: {""})
            update {bool} -- update isomesh or not? (default: {False})
        """
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
    CUBE = Cube()
    CUBE.load_header('test.cube')
    print(CUBE)
    CUBE.load_body()
    CUBE.make_color_voxel(update=True)
    CUBE.make_emission_voxel(update=True)
