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
ensemble = glob.glob("/home/asiera/bluebottle/cases/shear-laura/*")

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
nt = 50
t_end = 0.5

####################################################################
# store particle position data at initial time t_init
####################################################################

timeseries = bb.init(ensemble[-1] + "/output")
bb.open(timeseries[0])  # using time step 0 for now
(XX, YY, ZZ) = bb.read_part_position()
X_init = numpy.array(XX)
Y_init = numpy.array(YY)
Z_init = numpy.array(ZZ)
T_init = bb.read_time()
np = len(X_init)  # also store particle number
bb.close()

# create particle lists
# last particle position
X0 = numpy.zeros(np)
Y0 = numpy.zeros(np)
Z0 = numpy.zeros(np)
# current particle position
X = numpy.zeros(np)
Y = numpy.zeros(np)
Z = numpy.zeros(np)
# current particle velocity
U = numpy.zeros(np)
V = numpy.zeros(np)
W = numpy.zeros(np)
# mean square displacements
xMSD = numpy.zeros((nt,np))
yMSD = numpy.zeros((nt,np))
zMSD = numpy.zeros((nt,np))
# time difference (t - T_init)
T = numpy.zeros(nt)

##################################
# compute mean square displacement
##################################

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
    print("\rrealization", rcount, "of", len(ensemble), ": time =", time, "of",
      round(t_end, 2), "(" + str(round(percent)) + "%)    ", end='')
    sys.stdout.flush()

    # open and process each time step in each realization
    bb.open(time)
    t = bb.read_time()
    if t < t_end: # only read until reach end time of shortest simulation
      [X, Y, Z] = bb.read_part_position()
      [xMSD[tind][:], yMSD[tind][:], zMSD[tind][:]] = bb.msd(X, Y, Z,
        X_init, Y_init, Z_init)
      T[tind] = t
      tind = tind + 1

      # store current positions
      X0 = numpy.copy(X)
      Y0 = numpy.copy(Y)
      Z0 = numpy.copy(Z)

      bb.close()
    else:
      bb.close()
      break
    percent = percent + dpercent

# add an extra linebreak for better stdout
print()

# plot averaged velocities
plt.plot(T,bb.part_mean(xMSD.transpose()[:]),label='xMSD')
#plt.plot(T,V,label='V(t)')
#plt.plot(T,W,label='W(t)')
#plt.legend()
plt.xlabel('time')
plt.ylabel('mean square displacement')
plt.show()
