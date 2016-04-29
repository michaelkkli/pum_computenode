#!/usr/bin/python

import sys

print "len(sys.argv) is ", len(sys.argv)

if len(sys.argv)<3 :
	print "Usage: ", sys.argv[0], " input, output [, activate_scale_bar]"
	exit(1)

# Recorded script from Mayavi2
from numpy import array
try:
    engine = mayavi.engine
except NameError:
    from enthought.mayavi.api import Engine
    engine = Engine()
    engine.start()
if len(engine.scenes) == 0:
    engine.new_scene(size=(512,512))
    #engine.new_scene(size=(1000,1000))
# ------------------------------------------- 
#vtkxml_file_reader = engine.open(u'/home/mike/vilnius/prog/computenode/src/test/test.vtp')
vtkxml_file_reader = engine.open(sys.argv[1])

from enthought.mayavi.modules.surface import Surface
surface = Surface()
engine.add_module(surface, obj=None)
scene = engine.scenes[0]
scene.scene.parallel_projection = True
scene.scene.background = (0.49803921568627452, 0.49803921568627452, 0.49803921568627452)
camera_light = engine.scenes[0].scene.light_manager.lights[0]
camera_light.elevation = 0.0
camera_light.azimuth = 0.0
camera_light1 = engine.scenes[0].scene.light_manager.lights[1]
camera_light1.elevation = 0.0
camera_light1.azimuth = 0.0
camera_light1.intensity = 1.0
camera_light1.activate = False
camera_light2 = engine.scenes[0].scene.light_manager.lights[2]
camera_light2.elevation = 0.0
camera_light2.azimuth = 0.0
camera_light2.intensity = 1.0
camera_light2.activate = False
scene.scene.light_manager.light_mode = 'vtk'
if (len(sys.argv)>3) and (sys.argv[3]=="1"):
	module_manager = vtkxml_file_reader.children[0]
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
	module_manager.scalar_lut_manager.scalar_bar.reference_count = 4
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
	module_manager.scalar_lut_manager.show_scalar_bar = True
	module_manager.scalar_lut_manager.use_default_name = False
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
	module_manager.scalar_lut_manager.scalar_bar.title = 'Percentage'
	module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
	module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
	module_manager.scalar_lut_manager.data_range = array([   0.,  100.])
	module_manager.scalar_lut_manager.data_name = u'Percentage'

scene.scene.magnification = 1
scene.scene.save(sys.argv[2], size=(512,512))

#scene.scene.save(u'/home/mike/vilnius/prog/computenode/src/test/snapshot.png')
#module_manager.scalar_lut_manager.scalar_bar.position2 = array([ 0.8 ,  0.17])
#module_manager.scalar_lut_manager.scalar_bar.position = array([ 0.1 ,  0.01])
#module_manager.scalar_lut_manager.data_range = array([   0.,  100.])


## From website on Mayavi2 command line arguments.
#$ mayavi2 -d ParametricSurface -s "function='dini'" -m Surface \
#  -s "module_manager.scalar_lut_manager.show_scalar_bar = True" \
#  -s "scene.isometric_view()" -s "scene.save('snapshot.png')"
#
#$ mayavi2 -d heart.vtk -m Axes -m Outline -m GridPlane \
#  -m ContourGridPlane -m IsoSurface
#
#$ mayavi2 -d fire_ug.vtu -m Axes -m Outline -m VectorCutPlane \
#  -f MaskPoints -m Glyph
