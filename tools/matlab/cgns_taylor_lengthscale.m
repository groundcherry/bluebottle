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

function [L t] = cgns_taylor_lengthscale(casename, ts, te, loc)
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
xl = x(end,1,1) - x(1,1,1);
dx = x(2,1,1) - x(1,1,1);
yl = y(1,end,1) - y(1,1,1);
dy = y(1,2,1) - y(1,1,1);
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

lambda_gx = zeros(nt, 1);
lambda_gy = zeros(nt, 1);
t = zeros(nt, 1);

kx = (0:Nx-1)';
ky = (0:Ny-1)';

for i = 1:nt
  % get interpolated fields
  [u v w] = cgns_read_flow_vel_plane_z(casename, tstr{n(i)}, loc);

  %X(1,:) = -6+(12/64):12/64:6-(12/64);
  %u(1,:) = sin(4*2*pi*X/12);
  %figure(1)
  %plot(X,u)
  %xlabel('x')
  %ylabel('phi(u)')
  %
  %Y(:,1) = -6+(12/64):12/64:6-(12/64);
  %v(:,1) = sin(2*2*pi*Y/12);
  %figure(2)
  %plot(Y,v)
  %xlabel('y')
  %ylabel('phi(v)')
  %
  %phi(:,:) = v*u;
  %figure(3)
  %surf(phi)
  %xlabel('x')
  %ylabel('y')
  %zlabel('phi(x,y)')

  % do FFT
  % These need to be inverted
  U = fft(u',Nx,2)/Nx;
  V = fft(v',Ny,1)/Ny;

  % do convolution; note the factor of two (see Numerical Recipes III p. 603)
  ex = U.*conj(U);
  ey = V.*conj(V);

  % do inverse FFT
  Rx = ifft(ex,Nx,2)*Nx;
  Ry = ifft(ey,Ny,1)*Ny;

  % divide by mean square velocity for this time step
  umean = mean(mean(u.*u));
  vmean = mean(mean(v.*v));

  Rx = Rx / umean;
  Ry = Ry / vmean;

  Ex = mean(Rx,1)';
  Ey = mean(Ry,2);

  fpp2x = (Ex(3)-2*Ex(2)+Ex(1))/dx/dx;
  fpp1x = (Ex(2)-Ex(1))/dx/dx; % assume even function
  fpp0x = fpp1x+0.5*dx*(fpp2x-fpp1x)/dx;
  lambda_gx(i) = 1/sqrt(-0.5*fpp0x) / sqrt(2);
  fpp2y = (Ey(3)-2*Ey(2)+Ey(1))/dy/dy;
  fpp1y = (Ey(2)-Ey(1))/dy/dy; % assume even function
  fpp0y = fpp1y+0.5*dy*(fpp2y-fpp1y)/dy;
  lambda_gy(i) = 1/sqrt(-0.5*fpp0y) / sqrt(2);

  t(i) = tnum(i);
end

L = 0.5*(lambda_gx + lambda_gy);

%figure(4)
plot(t,L)
%xlabel('$k$','interpreter','latex')
%ylabel('$E(k)$','interpreter','latex')
