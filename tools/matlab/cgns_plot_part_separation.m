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

function cgns_plot_part_separation(casename, ts, te, n0, n1)
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
a = cgns_read_part_radius(casename,ts);

[np d] = size(x); % np = number of particles

% create data storage locations
X = zeros(np, nt);
Y = zeros(np, nt);
Z = zeros(np, nt);

% read time and particle position data
for i = 1:nt
  [X(:,i), Y(:,i), Z(:,i)] = cgns_read_part_position(casename, tstr{n(i)});
end

XX = X(n1+1,:) - X(n0+1,:);
for i = 1:nt
  if(XX(i) > (16 - a(n0) - a(n1))) XX(i)  = XX(i)  - 16;  %%%%%%%%%%%% use grid to get domain size
  elseif(XX(i)  < (-16 + a(n0) + a(n1))) XX(i)  = XX(i)  + 16;  %%%%%%%%%%%% use grid to get domain size
  end
end
YY = Y(n1+1,:) - Y(n0+1,:);
for i = 1:nt
  if(YY(i)  > (16 - a(n0) - a(n1))) YY(i)  = YY(i)  - 16;
  elseif(YY(i)  < (-16 + a(n0) + a(n1))) YY(i)  = YY(i)  + 16;  %%%%%%%%%%%% use grid to get domain size
  end
end
ZZ = Z(n1+1,:) - Z(n0+1,:);
for i = 1:nt
  if(ZZ(i)  > (16 - a(n0) - a(n1))) ZZ(i)  = ZZ(i)  - 16;
  elseif(ZZ(i)  < (-16 + a(n0) + a(n1))) ZZ(i)  = ZZ(i)  + 16;  %%%%%%%%%%%% use grid to get domain size
  end
end
dist = sqrt(XX.*XX+YY.*YY+ZZ.*ZZ)-a(n0)-a(n1);

% plot
figure
plot(tnum, dist)
xlabel('$t$', 'interpreter', 'latex')
ylabel('$h$', 'interpreter', 'latex')
