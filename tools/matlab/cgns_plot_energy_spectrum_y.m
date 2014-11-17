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

function [E k] = cgns_plot_energy_spectrum_y(casename, ts, te, loc)
% CGNS_PLOT_ENERGY_SPECTRUM_Y Plot the turbulent kinetic energy spectrum
%   computed at a plane normal to the y-axis at location loc.
%
%   [E k] = CGNS_PLOT_ENERGY_SPECTRUM_Y (CASENAME, TIME, LOC) plots the
%   turbulent kinetic energy computed at a plane normal to the y-axis at the
%   location LOC from the simulation CASENAME at time TIME.
%
%   Example:
%     cgns_plot_energy_spectrum_y('simulation', 3.14159, 4.2) will read the
%     appropriate output file located in 'simulation/output' at location
%     y = 4.2.

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
ys = y(1,1,1);
ye = y(1,end,1);

% check that loc is inside grid
if loc < ys | loc > ye
  disp(sprintf('loc must fall between %f and %f', ys, ye))
  return % quit
end

% get number of grid points for FFT
[sx sy sz] = size(y);
Nz = (sz - 1);
Nx = (sx - 1);
Ez = zeros(Nz, 1);
Ex = zeros(Nx, 1);

kz = (0:Nz-1)';
kx = (0:Nx-1)';

for i = 1:nt
  % get interpolated fields
  [u v w] = cgns_read_flow_vel_plane_y(casename, tstr{n(i)}, loc);

  % do FFT
  W = fft(w,Nz,2)/Nz;
  U = fft(u,Nx,1)/Nx;

% do convolution; note the factor of two (see Numerical Recipes III p. 603)
  ez = 2*W.*conj(W);
  ex = 2*U.*conj(U);

  % take off the mean velocity
  ez(:,1) = 0;
  ex(1,:) = 0;

  % accumulate the statistics
  Ez = Ez + mean(ez,1)';
  Ex = Ex + mean(ex,2);
end
kz = kz(1:length(kz)/2);
Ez = Ez(1:length(kz))/nt;
kx = kx(1:length(kx)/2);
Ex = Ex(1:length(kx))/nt;

E = 0.5 * (Ez + Ex);
k = kz;

loglog(k,E)
xlabel('$k$','interpreter','latex')
ylabel('$E(k)$','interpreter','latex')
