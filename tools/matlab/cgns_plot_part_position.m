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

function cgns_plot_part_position(casename, ts, te)
% CGNS_PLOT_PART_POSITION  Plot the particle positions from a BLUEBOTTLE-
%   generated CGNS file.
%
%   CGNS_PLOT_PART_POSITION(CASENAME, TS, TE) plots the particle positions
%   from the simulation CASENAME at all times TS <= T <= TE
%
%   Example:
%     cgns_plot_part_position('simulation', 0, 1) will plot the position
%     history for the CASENAME = 'simulation' for 0 <= T <= 1

% read part time
[tstr tnum] = cgns_read_part_time(casename);
n = 1:1:length(tnum);

% 'cut' the time array between the ts and te values, inclusive
i = find(tnum<ts | tnum>te);
n(i) = [];
tnum(i) = [];
ts = tstr{n(1)};


% calculate the number of timesteps
nt = length(n);

% find the number of particles using the initial timestep
[x, y, z] = cgns_read_part_position(casename, ts);

[np d] = size(x); % np = number of particles

% create data storage locations
X = zeros(np, nt);
Y = zeros(np, nt);
Z = zeros(np, nt);
R = zeros(np, nt);

% read time and particle position data
for i = 1:nt
  [X(:,i), Y(:,i), Z(:,i)] = cgns_read_part_position(casename, tstr{n(i)});
  R(:,i) = cgns_read_part_radius(casename, tstr{n(i)});
end

% plot
figure
plot(tnum, X)
xlabel('$t$', 'interpreter', 'latex')
ylabel('$x$', 'interpreter', 'latex')

figure
plot(tnum, Y)
xlabel('$t$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')

figure
plot(tnum, Z)
xlabel('$t$', 'interpreter', 'latex')
ylabel('$z$', 'interpreter', 'latex')

for i=1:np
  for j=1:np
    if i~=j
      xx = X(i,:)-X(j,:);
      yy = Y(i,:)-Y(j,:);
      zz = Z(i,:)-Z(j,:);
      dist(i+np*(j-1),:) = sqrt(xx.*xx+yy.*yy+zz.*zz) - R(i,:) - R(j,:);
    end
  end
end

figure
plot(tnum, dist)
xlabel('$t$', 'interpreter', 'latex')
ylabel('$r$', 'interpreter', 'latex')
