#!/usr/bin/env python

################################################################################
################################## BLUEBOTTLE ##################################
################################################################################
#
#  Copyright 2012 - 2016 Adam Sierakowski, The Johns Hopkins University
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#  Please contact the Johns Hopkins University to use Bluebottle for
#  commercial and/or for-profit applications.
################################################################################

import glob, sys
import bluebottle_particle_reader as bb
import numpy
import matplotlib.pyplot as plt

# get all base directories
ensemble = glob.glob("/scratch/users/asierak1@jhu.edu/shear-laura/*")
#ensemble = glob.glob("/scratch/users/asierak1@jhu.edu/shear-laura-test/*")

############################################################
# visit each realization to find minimum simulation duration
############################################################

# maximum time of shortest realization
t_end = 10e10;

# number of time outputs
nt = 0;

for realization in ensemble:
  timeseries = bb.init(realization + "/output")

  # open final output in timeseries and read time for comparison
  bb.open(timeseries[-1])
  t_tmp = bb.read_time()
  if t_tmp < t_end:
    t_end = t_tmp
    nt = len(timeseries)

  # close this output
  bb.close()

# overwrite number of time outputs to read (for testing)
#nt = 5000
#t_end = 50.0

####################################################################
# average particle velocity over all particles as a function of time
####################################################################

# open, read, and close first realization to find particle number
timeseries = bb.init(ensemble[-1] + "/output")
bb.open(timeseries[0])
np = len(bb.read_part_position()[0]) # [0] ==> x-component
bb.close()

# create particle lists
U = numpy.zeros(nt)
V = numpy.zeros(nt)
W = numpy.zeros(nt)
T = numpy.zeros(nt)

# realization counter
rcount = 0
percent = 0
dpercent = 100. / len(ensemble) / nt
for realization in ensemble:
  rcount = rcount + 1
  timeseries = bb.init(realization + "/output")

  tind = 0 # time index
  # read all time outputs in timeseries
  for time in timeseries:
    # tell user where we are
    #print("\rrealization", rcount, "of", len(ensemble), ": time =", time, "of",
    #  round(t_end, 2), "(" + str(round(percent)) + "%)    ", end='')
    print("realization", rcount, "of", len(ensemble), ": time =", time, "of",
      round(t_end, 2), "(" + str(round(percent)) + "%)    ")
    sys.stdout.flush()

    # open and process each time step in each realization
    f = bb.open(time)
    if f != None:
      t = bb.read_time()
      if t < t_end: # only read until reach end time of shortest simulation
        [u, v, w] = bb.read_part_velocity()
        U[tind] = U[tind] + bb.part_mean(u)
        V[tind] = V[tind] + bb.part_mean(v)
        W[tind] = W[tind] + bb.part_mean(w)
        T[tind] = t
        tind = tind + 1
        bb.close()
      else:
        bb.close()
        break
    percent = percent + dpercent

# add an extra linebreak for better stdout
print()

# finish averaging
U = U / len(ensemble)
V = V / len(ensemble)
W = W / len(ensemble)

# plot averaged velocities
plt.plot(T,U,label='U(t)')
plt.plot(T,V,label='V(t)')
plt.plot(T,W,label='W(t)')
plt.legend()
plt.xlabel('time')
plt.ylabel('particle velocity')
fig = plt.gcf()
fig.set_size_inches(10,7.5)
fig.savefig('mean_velocity.png')
