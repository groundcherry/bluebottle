#!/home/asiera/sandbox/ParaView-4.2.0-Linux-64bit/bin/pvpython

# README
#
# Dependencies (I have tested it with these; others may also work):
#  - ParaView-4.2.0
#  - ffmpeg release/2.8 (git://source.ffmpeg.org/ffmpeg.git)
#
# Usage:
#  - Change the path given in Line 1 to point to pvpython inside your Paraview
#    installation
#  - Use ParaView to set up the scene as you like
#  - Make sure to rename the filter containing the data as 'flow' and 'part'
#  - Save the ParaView state file (File->Save State)
#  - Create a directory for storing the image files
#  - Run animate-cgns.py and follow the instructions to generate the animation
#
# Don't be afraid to edit this script for your own purposes (e.g., if you don't
# need to read a 'part' file, comment out those lines below.
#
# 2015 Adam Sierakowski sierakowski@jhu.edu

import sys
from paraview.simple import *
import glob, os
import subprocess

# TODO accept starting and ending time command line input

print "CGNS animation utility"
print ""
root = raw_input("Simulation root: ")
if not root.endswith('/'):
  root = root + '/'
state = raw_input("ParaView state file (.pvsm): ")
if not state.endswith('.pvsm'):
  state = state + '.pvsm'
ts = raw_input("Animation start time: ")
te = raw_input("Animation end time: ")
img = raw_input("Image output directory: ")
if not img.endswith('/'):
  img = img + '/'
anim = raw_input("Animation output name (include extension): ")
fps = raw_input("Animation frame rate: ")

print "\nSummary:"
print "  Simulation root: " + root
print "  ParaView state file: " + root + state
print "  Animation time: " + ts + " to " + te
print "  Image output directory: " + root + img
print "  Animation output name: " + root + anim
print "  Animation frame rate: " + fps

print ""
run = raw_input("Continue? (y/N) ")
print ""
if  run == 'y' or run == 'Y':
  # load state file
  # make sure to build state file using the filter names 'flow' and 'part'
  LoadState(state)
  view = GetRenderView()
  flow = FindSource('flow')
  if flow == None:
    print "Make sure to build the state file using the filter name \'flow\'"
  part = FindSource('part')
  if part == None:
    print "Make sure to build the state file using the filter name \'part\'"

  view.WriteImage(img + "tmp.png", "vtkPNGWriter", 1)
  os.remove(img + "tmp.png")

  # go through all files
  for fname in sorted(glob.glob(root + "/output/flow*.cgns")):
    time = fname.split('/')[-1]

    if time.endswith('.cgns'):
      time = time[:-5]
    if time.startswith('flow-'):
      time = time[5:]

    if float(time) >= float(ts) and float(time) <= float(te):
      # change to file given by time
      flow.FileName = root + "/output/flow-" + time + ".cgns"
      flow.FileNameChanged()
      part.FileName = root + "/output/part-" + time + ".cgns"
      part.FileNameChanged()

      print "Saving image for t = " + time
      
      # save screen shot
      view.WriteImage(img + "img-" + time + ".png", "vtkPNGWriter", 1)

  # stitch together using ffmpeg
  print ""
  print "Stitching images together into " + anim
  print ""
  pts = "setpts=" + str(30./float(fps)) + "*PTS"
  subprocess.call(["ffmpeg", "-r", "30", "-f", "image2", "-pattern_type",
    "glob", "-i", img + "*.png", "-qscale:v", "4", "-vcodec",
    "msmpeg4v2", "-r", "30", "-vf", pts, anim])
