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

function [x, y, z] = cgns_read_part_position(casename, time)
% CGNS_READ_PART_POSITION Read the particle positions from a BLUEBOTTLE-
%   generated CGNS file.
%
%   [x, y, z] = CGNS_READ_PART_POSITION(CASENAME, TIME) reads the particle
%   positions from the simulation CASENAME at time TIME. Each of the
%   position components is an array representing all of the particles in
%   the simulation.
%
%   Example:
%     cgns_read_part_position('simulation', 3.14159) will read the
%     appropriate output file located in 'simulation/output

% determine input type of 'time'

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

xsol = '/Base/Zone0/GridCoordinates/CoordinateX/ data';
ysol = '/Base/Zone0/GridCoordinates/CoordinateY/ data';
zsol = '/Base/Zone0/GridCoordinates/CoordinateZ/ data';

x = h5read(path, xsol);
y = h5read(path, ysol);
z = h5read(path, zsol);

