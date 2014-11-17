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

function [x, y, z] = cgns_read_grid(casename)
% CGNS_READ_GRID  Read the discretization grid from a BLUEBOTTLE-generated CGNS
%   file.
%
%   [x, y, z] = CGNS_READ_GRID(CASENAME) reads the flow velocity field from the
%   simulation CASENAME. Each component holds a three-dimensional array
%   whose values provide the value of the corresponding Cartesian vector at
%   each position. For example, the coordinates of a single point in 3D space
%   can be constructed by: x(i,j,k)*i_hat + y(i,j,k)*j_hat + z(i,j,k)*k_hat.
%
%   Example:
%     cgns_read_grid('simulation') will read the appropriate output file located
%     in 'simulation/output'.

path = casename;
path = [path '/output/grid.cgns'];

xsol = '/Base/Zone0/GridCoordinates/CoordinateX/ data';
ysol = '/Base/Zone0/GridCoordinates/CoordinateY/ data';
zsol = '/Base/Zone0/GridCoordinates/CoordinateZ/ data';

x = h5read(path, xsol);
y = h5read(path, ysol);
z = h5read(path, zsol);
