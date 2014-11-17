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

function [E k] = cgns_plot_energy_spectrum_x(casename, ts, te, loc)
% CGNS_PLOT_ENERGY_SPECTRUM_X Plot the turbulent kinetic energy spectrum
%   computed at a plane normal to the x-axis at location loc.
%
%   [E k] = CGNS_PLOT_ENERGY_SPECTRUM_X (CASENAME, TIME, LOC) plots the
%   turbulent kinetic energy computed at a plane normal to the x-axis at the
%   location LOC from the simulation CASENAME at time TIME.
%
%   Example:
%     cgns_plot_energy_spectrum_x('simulation', 3.14159, 4.2) will read the
%     appropriate output file located in 'simulation/output' at location
%     x = 4.2.

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
xs = x(1,1,1);
xe = x(end,1,1);

% check that loc is inside grid
if loc < xs | loc > xe
  disp(sprintf('loc must fall between %f and %f', xs, xe))
  return % quit
end

% get number of grid points for FFT
[sx sy sz] = size(x);
Ny = (sy - 1);
Nz = (sz - 1);
Ey = zeros(Ny, 1);
Ez = zeros(Nz, 1);

ky = (0:Ny-1)';
kz = (0:Nz-1)';

for i = 1:nt
  % get interpolated fields
  [u v w] = cgns_read_flow_vel_plane_x(casename, tstr{n(i)}, loc);

  % do FFT
  V = fft(v,Ny,2)/Ny;
  W = fft(w,Nz,1)/Nz;

  % do convolution; note the factor of two (see Numerical Recipes III p. 603)
  ey = 2*V.*conj(V);
  ez = 2*W.*conj(W);

  % take off the mean velocity
  ey(:,1) = 0;
  ez(1,:) = 0;

  % accumulate the statistics
  Ey = Ey + mean(ey,1)';
  Ez = Ez + mean(ez,2);
end
ky = ky(1:length(ky)/2);
Ey = Ey(1:length(ky))/nt;
kz = kz(1:length(kz)/2);
Ez = Ez(1:length(kz))/nt;

E = 0.5 * (Ey + Ez);
k = ky;

loglog(k,E)
xlabel('$k$','interpreter','latex')
ylabel('$E(k)$','interpreter','latex')
