#!/usr/bin/env python

# bluebottle_particle_reader python module example code

import sys, getopt
import numpy as np
import bluebottle_particle_reader as bb

# initialize the reader
times = bb.init("/home/asiera/bluebottle/sim/output")

# visit all outputted time values
for time in times:
  # open the CGNS file for this particular output time
  bb.open(time)

  # read the CGNS file
  t = bb.read_time()
  (x,y,z) = bb.read_part_position()
  (u,v,w) = bb.read_part_velocity()

  print("t =", t)

  # close the CGNS file
  bb.close()
