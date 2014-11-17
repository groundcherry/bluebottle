%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLUEBOTTLE-1.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Copyright 2012 - 2014 Adam Sierakowski, The Johns Hopkins University
% 
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
% 
%       http://www.apache.org/licenses/LICENSE-2.0
% 
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.
% 
%   Please contact the Johns Hopkins University to use Bluebottle for
%   commercial and/or for-profit applications.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [udot, vdot, wdot] = cgns_read_part_acc(casename, time)
% CGNS_READ_PART_ACC Read the particle accelerations from a BLUEBOTTLE-
%   generated CGNS file.
%
%   [udot, vdot, wdot] = CGNS_READ_PART_ACC(CASENAME, TIME) reads the particle
%   accelerations from the simulation CASENAME at time TIME. Each of the
%   acceleration components is an array representing all of the particles
%   in the simulation.
%
%   Example:
%     cgns_read_part_acc('simulation', 3.14159) will read the appropriate
%     output file located in 'simulation/output'

if isa(time, 'double') == 1
  % find the directory contents
  path = [casename '/output'];
  od = cd(path);
  contents = dir;
  % Get names of files, remove first two (.., .), find flows/parts
  names = {contents.name}';
  names(1:2) = [];
  check = find(strncmp(names, 'part', 4));
  % Get sigfigs
  sigfigs = names{check(1)};
  sigfigs = length(sigfigs(6:end-5)) - 2;
  t_format = sprintf('%%1.%df', sigfigs);
  tt = sprintf(t_format, time);
  cd(od);
elseif isa(time, 'char') == 1
    tt = time;
end

path = [casename '/output/part-' tt '.cgns'];

usol = '/Base/Zone0/Solution/AccelerationX/ data';
vsol = '/Base/Zone0/Solution/AccelerationY/ data';
wsol = '/Base/Zone0/Solution/AccelerationZ/ data';

udot = h5read(path, usol);
vdot = h5read(path, vsol);
wdot = h5read(path, wsol);
