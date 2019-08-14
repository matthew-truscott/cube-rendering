import bpy
import molecule
import json
import utils_blender as ub
import numpy as np

class CPKData():
    def __init__(self):
        self.colors = {
            1: (4, 4, 4),
            6: (0.5, 0.5, 0.5),
            7: (2, 2, 3.56),
            8: (4, 0, 0)
        }

        with open('ions.json', 'r') as f:
            self.radii = json.load(f)

        self.atom_scale = 0.24
        self.bond_width = 0.24

def edit_atom_material(molecule, diffuse=0.9, specular=0.0, hardness=50):
    # in the future, add routines to selectively add different materials to atoms
    for obj in molecule.rendered:
        main_material = obj.data.materials[0]
        main_material.diffuse_intensity = diffuse
        main_material.specular_intensity = specular
        main_material.specular_hardness = hardness

def draw_molecule(cube, bonds=False):
    molecule = cube.molecule

    # fix positions if the field was rolled
    if cube.settings.roll:
        cube.molecule.transform(-1.0 * cube.field.transform)

    cpkdata = CPKData()
    atoms = []
    # draw spheres
    print('drawing atoms...')
    for i, p in enumerate(molecule.m_positions):
        bpy.ops.mesh.primitive_uv_sphere_add(location=(p[0], p[1], p[2]), segments=64, ring_count=32)
        if molecule.a_species[i] in cpkdata.colors:
            atom_color = cpkdata.colors[molecule.a_species[i]]
        else:
            atom_color = (0.8, 0.4, 0.4)
        smat = ub.simpleMaterial(atom_color)
        sobj = bpy.context.active_object
        sobj.data.materials.append(smat)
        # scale accordingly
        uff = float(cpkdata.radii['{}'.format(molecule.a_species[i])]) * cpkdata.atom_scale
        sobj.scale = (uff, uff, uff)
        atoms.append(sobj)

    molecule.add_rendered(atoms)

    if not bonds:
        return molecule

    if not molecule.bondlist:
        molecule.create_bonds()

    if not molecule.bondlist:
        return molecule

    # draw bonds, split each into two cylinders to match colors and stuff
    print('drawing bonds...')
    for bond in molecule.bondlist:
        # first find the midway point, the 'position' of the cylinder
        pointA = molecule.m_positions[bond[0]]
        pointB = molecule.m_positions[bond[1]]
        pointM = (pointA + pointB) / 2.0
        pointO = [0, 0, 0]
        vecA = [0, 0, 0]

        for i in range(3):
            pointO[i] = (pointA[i] + pointM[i]) / 2.0
            vecA[i] = (pointA[i] - pointM[i])

        # next, get rotation matrix
        vecA = np.array(vecA)
        normA = np.linalg.norm(vecA)
        vecA = vecA / normA
        vecB = np.array([0, 0, 1])
        vecC = np.cross(vecA, vecB)
        dot = np.dot(vecA, vecB)
        skew = np.array([[0, -vecC[2], vecC[1]],
                         [vecC[2], 0, -vecC[0]],
                         [-vecC[1], vecC[0], 0]])
        rotator = np.identity(3) + skew + (np.matmul(skew, skew) * (1.0 / (1.0 + dot)))

        # convert to 4-matrix
        rot4 = np.zeros((4, 4), dtype=float)
        rot4[0:3, 0:3] = rotator
        rot4[3, 3] = 1

        # now use this rotation matrix to rotate the bond
        bpy.ops.mesh.primitive_cylinder_add(location=(pointO[0], pointO[1], pointO[2]))
        sobj1 = bpy.context.active_object

        # scale accordingly
        scale4 = np.eye(4, dtype=float)
        scale4[0, 0] = cpkdata.bond_width
        scale4[1, 1] = cpkdata.bond_width
        scale4[2, 2] = normA / 2

        # apply these transformations
        sobj1.data.transform(scale4)
        sobj1.data.transform(rot4)
        sobj1.data.update()

        # color accordingly
        atom1 = molecule.a_species[bond[0]]
        color1 = cpkdata.colors[atom1]
        smat1 = ub.simpleMaterial(color1)
        sobj1.data.materials.append(smat1)

        # second cylinder is exactly the same as above (will skip comments)
        for i in range(3):
            pointO[i] = (pointM[i] + pointB[i]) / 2.0
            vecA[i] = (pointM[i] - pointB[i])
        vecA = np.array(vecA)
        normA = np.linalg.norm(vecA)
        vecA = vecA / normA
        vecB = np.array([0, 0, 1])
        vecC = np.cross(vecA, vecB)
        dot = np.dot(vecA, vecB)
        skew = np.array([[0, -vecC[2], vecC[1]],
                         [vecC[2], 0, -vecC[0]],
                         [-vecC[1], vecC[0], 0]])
        rotator = np.identity(3) + skew + (np.matmul(skew, skew) * (1.0 / (1.0 + dot)))
        rot4 = np.zeros((4, 4), dtype=float)
        rot4[0:3, 0:3] = rotator
        rot4[3, 3] = 1
        bpy.ops.mesh.primitive_cylinder_add(location=(pointO[0], pointO[1], pointO[2]))
        sobj2 = bpy.context.active_object
        scale4 = np.eye(4, dtype=float)
        scale4[0, 0] = cpkdata.bond_width
        scale4[1, 1] = cpkdata.bond_width
        scale4[2, 2] = normA / 2
        sobj2.data.transform(scale4)
        sobj2.data.transform(rot4)
        sobj2.data.update()
        atom2 = molecule.a_species[bond[1]]
        color2 = cpkdata.colors[atom2]
        smat2 = ub.simpleMaterial(color2)
        sobj2.data.materials.append(smat2)

    return molecule

