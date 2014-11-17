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

function [u, v, w, t] = cgns_read_point_vel(casename, ts, te, x, y, z)
% CGNS_READ_POINT_VEL  Read the flow velocity time series from a BLUEBOTTLE-
%   generated CGNS file at a point.
%
%   CGNS_READ_POINT_VEL(CASENAME, TS, TE, X, Y, Z) plots the particle velocity
%   time series at point (X, Y, Z) from the simulation CASENAME between TS
%   and TE.
%
%   Example:
%     cgns_read_point_vel('simulation',0,1,2,3,4) will plot the velocity time
%     series for the CASENAME = 'simulation' between t = 0 and 1 at
%     point (2,3,4).

% read part time
[tstr tnum] = cgns_read_flow_time(casename);
n = 1:1:length(tnum);

% 'cut' the time array between the ts and te values, inclusive
i = find(tnum<ts | tnum>te);
i(1) = []; % include the te value
n(i) = [];
tnum(i) = [];
ts = tstr{n(1)};

% calculate the number of timesteps
nt = length(n);

% create data storage locations
t = tnum;
u = zeros(nt, 1);
v = zeros(nt, 1);
w = zeros(nt, 1);

% get grid
[X Y Z] = cgns_read_grid(casename);
[xn yn zn] = size(X);
xn = xn - 1;
yn = yn - 1;
zn = zn - 1;
xs = X(1,1,1);
xe = X(end,1,1);
ys = Y(1,1,1);
ye = Y(1,end,1);
zs = Z(1,1,1);
ze = Z(1,1,end);
xl = xe - xs;
yl = ye - ys;
zl = ze - zs;
dx = xl / xn;
dy = yl / yn;
dz = zl / zn;

% set up interpolation
i = floor((x - xs)/dx + 0.5);
j = floor((y - ys)/dy + 0.5);
k = floor((z - zs)/dz + 0.5);

% read particle position data
for l = 1:nt
  [U, V, W] = cgns_read_flow_vel(casename, tstr{n(l)});
  % do interpolation
  if i == 0 % use one-sided extrapolation
    i = 1;
  elseif i == xn % use one-sided extrapolation
    i = xn-1;
  end
  if j == 0 % use one-sided extrapolation
    j = 1;
  elseif j == yn % use one-sided extrapolation
    j = yn-1;
  end
  if k == 0 % use one-sided extrapolation
    k = 1;
  elseif k == zn % use one-sided extrapolation
    i = kn-1;
  end
  % otherwise, use central difference
  xx = (i-0.5)*dx+xs;
  yy = (j-0.5)*dy+ys;
  zz = (k-0.5)*dz+zs;
  u00 = U(i,j,k) + (x-xx)/dx * (U(i+1,j,k)-U(i,j,k));
  u01 = U(i,j+1,k) + (x-xx)/dx * (U(i+1,j+1,k)-U(i,j+1,k));
  u10 = U(i,j,k+1) + (x-xx)/dx * (U(i+1,j,k+1)-U(i,j,k+1));
  u11 = U(i,j+1,k+1) + (x-xx)/dx * (U(i+1,j+1,k+1)-U(i,j+1,k+1));
  u0 = u00 + (y-yy)/dy * (u01-u00);
  u1 = u10 + (y-yy)/dy * (u11-u10);
  u(l) = u0 + (z-zz)/dz * (u1-u0);
  v00 = V(i,j,k) + (x-xx)/dx * (V(i+1,j,k)-V(i,j,k));
  v01 = V(i,j+1,k) + (x-xx)/dx * (V(i+1,j+1,k)-V(i,j+1,k));
  v10 = V(i,j,k+1) + (x-xx)/dx * (V(i+1,j,k+1)-V(i,j,k+1));
  v11 = V(i,j+1,k+1) + (x-xx)/dx * (V(i+1,j+1,k+1)-V(i,j+1,k+1));
  v0 = v00 + (y-yy)/dy * (v01-v00);
  v1 = v10 + (y-yy)/dy * (v11-v10);
  v(l) = v0 + (z-zz)/dz * (v1-v0);
  w00 = W(i,j,k) + (x-xx)/dx * (W(i+1,j,k)-W(i,j,k));
  w01 = W(i,j+1,k) + (x-xx)/dx * (W(i+1,j+1,k)-W(i,j+1,k));
  w10 = W(i,j,k+1) + (x-xx)/dx * (W(i+1,j,k+1)-W(i,j,k+1));
  w11 = W(i,j+1,k+1) + (x-xx)/dx * (W(i+1,j+1,k+1)-W(i,j+1,k+1));
  w0 = w00 + (y-yy)/dy * (w01-w00);
  w1 = w10 + (y-yy)/dy * (w11-w10);
  w(l) = w0 + (z-zz)/dz * (w1-w0);
end
