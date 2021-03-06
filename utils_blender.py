"""UTILS_BLENDER MODULE

Adapted from https://github.com/njanakiev/blender-scripting

"""
import colorsys
import os

import bpy
import bmesh
from math import sin, cos, pi
tau = 2.0 * pi


def set_background(color):
    if 'World' in bpy.data.worlds:
        breakpoint()
        bpy.data.worlds['World'].horizon_color = color


def removeObject(obj):
    if obj.type == 'MESH':
        if obj.data.name in bpy.data.meshes:
            bpy.data.meshes.remove(obj.data)
        if obj.name in bpy.context.scene.objects:
            bpy.context.scene.objects.unlink(obj)
        bpy.data.objects.remove(obj)
    else:
        raise NotImplementedError('Other types not implemented yet besides \'MESH\'')


def trackToConstraint(obj, target):
    constraint = obj.constraints.new('TRACK_TO')
    constraint.target = target
    constraint.track_axis = 'TRACK_NEGATIVE_Z'
    constraint.up_axis = 'UP_Y'

    return constraint


def target(origin=(0,0,0)):
    bpy.ops.object.empty_add(type="PLAIN_AXES", location=origin)
    obj = bpy.context.active_object
    return obj


def camera(origin, target=None, lens=45, clip_start=0.1, clip_end=200, type='PERSP', ortho_scale=6):
    bpy.ops.object.camera_add()
    camera = bpy.context.active_object
    camera.data.lens = lens
    camera.data.clip_start = clip_start
    camera.data.clip_end = clip_end
    camera.data.type = type
    camera.location = origin
    if type == 'ORTHO':
        camera.data.ortho_scale = ortho_scale

    if target: 
        trackToConstraint(camera, target)
    return camera

def get_scene():
    return bpy.context.scene


def lamp(origin, type='POINT', energy=1, color=(1,1,1), target=None):
    print('createLamp called')
    bpy.ops.object.light_add(type=type, location=origin)
    obj = bpy.context.active_object
    obj.data.energy = energy
    obj.data.color = color

    if target: 
        trackToConstraint(obj, target)
    return obj


def simpleScene(targetCoord, cameraCoord, sunCoord, lens=35):
    print('createSimpleScene called')

    tar = target(targetCoord)
    cam = camera(cameraCoord, tar, lens)
    sun = lamp(sunCoord, 'SUN', target=tar)

    return tar, cam, sun


def setAmbientOcclusion(ambient_occlusion=True, samples=5, blend_type='ADD'):
    # blend_type options: 'ADD', 'MULTIPLY'
    bpy.context.scene.world.light_settings.use_ambient_occlusion = ambient_occlusion
    bpy.context.scene.world.light_settings.ao_blend_type = blend_type
    bpy.context.scene.world.light_settings.samples = samples


def setSmooth(obj, level=None, smooth=True):
    if level:
        # Add subsurf modifier
        modifier = obj.modifiers.new('Subsurf', 'SUBSURF')
        modifier.levels = level
        modifier.render_levels = level

    # smooth surface
    mesh = obj.data
    for p in mesh.polygons:
        p.use_smooth = smooth


def cube(origin, size):
    bpy.ops.mesh.primitive_cube_add(location=origin)
    obj = bpy.context.active_object
    obj.scale(size, size, size)

    return obj


def rainbowLights(r=5, n=100, freq=2, energy=0.1):
    for i in range(n):
        t = float(i) / float(n)
        pos = (r * sin(tau * t), r * cos(tau * t), t* sin(freq * tau * t))

        # Create lamp
        bpy.ops.object.add(type='LAMP', location=pos)
        obj = bpy.context.object
        obj.data.type = 'POINT'

        # Apply gamma correction for Blender
        color = tuple(pow(c, 2.2) for c in colorsys.hsv_to_rgb(t, 0.6, 1))

        # Set HSV color and lamp energy
        obj.data.color = color
        obj.data.energy = energy


def removeAll(type=None):
    # possible type: 'MESH', 'CURVE', 'SURFACE', 'META', 'FONT', 'ARMATURE', 'LATTICE', 'EMPTY', 'CAMERA', 'LAMP'
    if type:
        bpy.ops.object.select_all(action='DESELECT')
        bpy.ops.object.select_by_type(type=type)
        bpy.ops.object.delete()
    # broken in blender 2.80
    else:
        override = bpy.context.copy()
        override['selected_objects'] = list(bpy.context.scene.objects)
        bpy.ops.object.delete(override)

def simpleMaterial(diffuse_color):
    mat = bpy.data.materials.new('Material')
    # Diffuse
    mat.diffuse_color = diffuse_color

    # Specular
    mat.specular_intensity = 0

    return mat


def volumeMaterial():
    mat = bpy.data.materials.new('VolumeMaterial')

    mat.type = 'VOLUME'
    mat.volume.scattering = 0.05
    mat.volume.reflection = 0.0


def falloffMaterial(diffuse_color):
    mat = bpy.data.materials.new('FalloffMaterial')

    # Diffuse
    mat.diffuse_shader = 'LAMBERT'
    mat.use_diffuse_ramp = True
    mat.diffuse_ramp_input = 'NORMAL'
    mat.diffuse_ramp_blend = 'ADD'
    mat.diffuse_ramp.elements[0].color = (1, 1, 1, 1)
    mat.diffuse_ramp_elements[1].color = (1, 1, 1, 0)
    mat.diffuse_color = diffuse_color
    mat.diffuse_intensity = 1.0

    # Specular
    mat.specular_intensity = 0.0

    # Shading
    mat.emit = 0.05
    mat.translucency = 0.2

    return mat


def voxel_texture(mat, filepath, field):
    tex = bpy.data.textures.new('VoxelTexture', 'VOXEL_DATA')
    tex.voxel_data.file_format = 'BLENDER_VOXEL'
    tex.voxel_data.filepath = filepath
    slot = mat.texture_slots.add()
    slot.texture = tex
    slot.texture_coords = 'ORCO'
    if field == 'color':
        slot.use_map_color_emission = True
    elif field == 'map':
        slot.use_map_emission = True

    return tex, slot


def colorRGB_256(color):
    return tuple(pow(float(c)/255.0, 2.2) for c in color)


def renderToFolder(renderFolder='rendering', renderName='render', resX=800, resY=800, resPercentage=100,
                   animation=False, frame_end=None):
    print('renderToFolder called')
    scn = bpy.context.scene
    scn.render.resolution_x = resX
    scn.render.resolution_y = resY
    scn.render.resolution_percentage = resPercentage
    if frame_end:
        scn.frame_end = frame_end

    print(bpy.context.space_data)

    # Check if script is executed inside Blender
    if bpy.context.space_data is None:
        # Specify folder to save rendering and check if it exists
        render_folder = os.path.join(os.getcwd(), renderFolder)
        if(not os.path.exists(render_folder)):
            os.mkdir(render_folder)

        if animation:
            # Render animation
            scn.render.filepath = os.path.join(render_folder, renderName)
            bpy.ops.render.render(animation=True)
        else:
            # Render still frame
            scn.render.filepath = os.path.join(render_folder, renderName + '.png')
            bpy.ops.render.render(write_still=True)


def save(filepath='untitled.blend'):
    bpy.ops.wm.save_as_mainfile(filepath=filepath, relative_remap=False)


def bmeshToObject(bm, name='Object'):
    mesh = bpy.data.meshes.new(name+'Mesh')
    bm.to_mesh(mesh)
    bm.free

    obj = bpy.data.objects.new(name, mesh)
    bpy.context.scene.objects.link(obj)
    bpy.context.scene.update()

    return obj
