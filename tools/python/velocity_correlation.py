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
#ensemble = glob.glob("/home/asiera/bluebottle/cases/shear-laura/*")
#ensemble = glob.glob("/home/asiera/bluebottle/tools/python/ensemble/*")
ensemble = glob.glob("/scratch/users/asierak1@jhu.edu/shear-laura/*")

# input some domain info (TODO: automate this)
Lx = 12
Ly = 12
Lz = 12
a = 1

timestart = 200
DT_out = 0.1

############################################################
# visit each realization to find minimum simulation duration
############################################################

# maximum time of shortest realization
t_end = 10e10;

# number of time outputs
nt = 0;

for realization in ensemble:
  timeseries = bb.init(realization + "/output")[int(timestart/DT_out):]

  # open final output in timeseries and read time for comparison
  bb.open(timeseries[-1])
  t_tmp = bb.read_time()
  if t_tmp < t_end:
    t_end = t_tmp
    nt = round(len(timeseries)/10) + 1

  # close this output
  bb.close()

# overwrite number of time outputs to read (for testing)
#nt = 10
#t_end = 21

####################################################################
# store particle position data at initial time t_init
####################################################################

timeseries = bb.init(ensemble[-1] + "/output")[int(timestart/DT_out):]
bb.open(timeseries[0])  # using time step 0 for now
(U_init, V_init, W_init) = bb.read_part_velocity()
T_init = bb.read_time()
np = len(U_init)  # also store particle number
bb.close()

# create particle lists
# current particle velocity
U = numpy.zeros(np)
V = numpy.zeros(np)
W = numpy.zeros(np)
# mean square displacements
uCOR = numpy.zeros((nt,np))
vCOR = numpy.zeros((nt,np))
wCOR = numpy.zeros((nt,np))
# mean for each particle
umean = numpy.zeros(np)
vmean = numpy.zeros(np)
wmean = numpy.zeros(np)
# counts
sumcount = numpy.zeros(nt)
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
  timeseries = bb.init(realization + "/output")[int(timestart/DT_out):]

  # read all time outputs in timeseries
  for time in timeseries[::10]:
    # tell user where we are
    print("realization", rcount, "of", len(ensemble), ": time =", time, "of",
      round(t_end, 2), "(" + str(round(percent)) + "%)    ")
    sys.stdout.flush()

    # set initial time
    bb.open(time)
    (U_init, V_init, W_init) = bb.read_part_velocity()
    T_init = bb.read_time()
    bb.close()
    if T_init < t_end:
      for time2 in timeseries[::10]:
        # open and process each time step in each realization
        bb.open(time2)
        t = bb.read_time()
        if t < t_end: # only read until reach end time of shortest simulation
          #dt = abs(t - T_init)
          dt = t - T_init
          if dt > 0:
            [U, V, W] = bb.read_part_velocity()

            [ucor, vcor, wcor] = bb.cor(U, V, W, U_init, V_init, W_init)

            dt_ind = int(round(dt / DT_out))
            print("\r  time =", time, " time2 =", time2, "dt =", dt, end='')

            uCOR[dt_ind] = uCOR[dt_ind] + ucor
            vCOR[dt_ind] = vCOR[dt_ind] + vcor
            wCOR[dt_ind] = wCOR[dt_ind] + wcor

            sumcount[dt_ind] = sumcount[dt_ind] + 1

            T[dt_ind] = dt

          bb.close()
        else:
          bb.close()
          break
      print()
    else:
      break
    percent = percent + dpercent

sumcount[sumcount == 0] = 1
# finish computing average
uCOR = numpy.divide(uCOR.transpose(),sumcount)
vCOR = numpy.divide(vCOR.transpose(),sumcount)
wCOR = numpy.divide(wCOR.transpose(),sumcount)

# add an extra linebreak for better stdout
print()

# plot averaged velocities
#plt.plot(T[:-1],bb.part_mean(xMSD.transpose()[:])[:-1],label='xMSD')
plt.plot(T[:-1],bb.part_mean(vCOR[:])[:-1],label='yMSD')
#plt.plot(T[:-1],bb.part_mean(zMSD.transpose()[:])[:-1],label='zMSD')
plt.legend()
plt.xlabel('time')
plt.ylabel('velocity correlation')
fig = plt.gcf()
fig.set_size_inches(10,7.5)
fig.savefig('velocity_correlation.png')
