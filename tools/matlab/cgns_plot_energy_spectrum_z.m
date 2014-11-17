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

function [E k] = cgns_plot_energy_spectrum_z(casename, ts, te, loc)
% CGNS_PLOT_ENERGY_SPECTRUM_Z Plot the turbulent kinetic energy spectrum
%   computed at a plane normal to the z-axis at location loc.
%
%   [E k] = CGNS_PLOT_ENERGY_SPECTRUM_Z (CASENAME, TIME, LOC) plots the
%   turbulent kinetic energy computed at a plane normal to the z-axis at the
%   location LOC from the simulation CASENAME at time TIME.
%
%   Example:
%     cgns_plot_energy_spectrum_z('simulation', 3.14159, 4.2) will read the
%     appropriate output file located in 'simulation/output' at location
%     z = 4.2.

% read flow time
[tstr tnum] = cgns_read_flow_time(casename);
n = 1:1:length(tnum);

% 'cut' the time array between the ts and te values, inclusive
i = find(tnum<ts | tnum>te);
n(i) = [];
tnum(i) = [];

% calculate the number of timesteps
nt = length(n);

% get grid for extents
[x y z] = cgns_read_grid(casename);
zs = z(1,1,1);
ze = z(1,1,end);

% check that loc is inside grid
if loc < zs | loc > ze
  disp(sprintf('loc must fall between %f and %f', zs, ze))
  return % quit
end

% get number of grid points for FFT
[sx sy sz] = size(z);
Nx = (sx - 1);
Ny = (sy - 1);
Ex = zeros(Nx, 1);
Ey = zeros(Ny, 1);

kx = (0:Nx-1)';
ky = (0:Ny-1)';

for i = 1:nt
  % get interpolated fields
  [u v w] = cgns_read_flow_vel_plane_z(casename, tstr{n(i)}, loc);

  % do FFT
  U = fft(u,Nx,2)/Nx;
  V = fft(v,Ny,1)/Ny;

  % do convolution; note the factor of two (see Numerical Recipes III p. 603)
  ex = 2*U.*conj(U);
  ey = 2*V.*conj(V);

  % take off the mean velocity
  ex(:,1) = 0;
  ey(1,:) = 0;

  % accumulate the statistics
  Ex = Ex + mean(ex,1)';
  Ey = Ey + mean(ey,2);
end

kx = kx(1:length(kx)/2);
Ex = Ex(1:length(kx))/nt;
ky = ky(1:length(ky)/2);
Ey = Ey(1:length(ky))/nt;

E = 0.5 * (Ex + Ey);
k = kx;

loglog(k,E)
xlabel('$k$','interpreter','latex')
ylabel('$E(k)$','interpreter','latex')
