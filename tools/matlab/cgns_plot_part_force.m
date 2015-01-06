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

function cgns_plot_part_force(casename, ts, te)
% CGNS_PLOT_PART_FORCE  Plot the particle forces from a BLUEBOTTLE-generated
%   CGNS file.
%
%   CGNS_PLOT_PART_FORCE(CASENAME, TS, TE, DT) plots the particle forces from
%   the simulation CASENAME at all times TS <= T <= TE
% 
%   If TS and TE do not exactly match the flow or part time given in the file
%   name, then the closest time step will be chosen
%
%   Example:
%     cgns_plot_part_force('simulation', 0.0025, 1) will plot the force history
%     for the CASENAME = 'simulation' for 0 <= T <= 1. If the real times 
%     are [0.002 0.0027..... 0.9958, 0.9978, 1.0001] then the force will
%     be plotted for 0.0027 through 0.9978

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
[fx, fy, fz] = cgns_read_part_force_total(casename, ts);

[np d] = size(fx); % np = number of particles

% create data storage locations
FX = zeros(np, nt);
FY = zeros(np, nt);
FZ = zeros(np, nt);
FXi = zeros(np, nt);
FYi = zeros(np, nt);
FZi = zeros(np, nt);
FXh = zeros(np, nt);
FYh = zeros(np, nt);
FZh = zeros(np, nt);

% read time and particle force data
for i = 1:nt
  [FX(:,i), FY(:,i), FZ(:,i)] = cgns_read_part_force_total(casename, tstr{n(i)});
  [FXi(:,i), FYi(:,i), FZi(:,i)] = cgns_read_part_force_interaction(casename, tstr{n(i)});
  [FXh(:,i), FYh(:,i), FZh(:,i)] = cgns_read_part_force_hydro(casename, tstr{n(i)});
end

% plot
figure
%plot(tnum, FX, tnum, FXi, '--', tnum, FXh, '.')
plot(tnum, FX, tnum, mean(FX,1), '-g.')%, tnum, FXi, '--', tnum, FXh, '.')
xlabel('$t$', 'interpreter', 'latex')
ylabel('$F_x$', 'interpreter', 'latex')
figure
%plot(tnum, FY, tnum, FYi, '--', tnum, FYh, '.')
plot(tnum, FY, tnum, mean(FY,1), '-g.')%, tnum, FYi, '--', tnum, FYh, '.')
xlabel('$t$', 'interpreter', 'latex')
ylabel('$F_y$', 'interpreter', 'latex')
figure
%plot(tnum, FZ, tnum, FZi, '--', tnum, FZh, '.')
plot(tnum, FZ, tnum, mean(FZ,1), '-g.')%, tnum, FZi, '--', tnum, FZh, '.')
xlabel('$t$', 'interpreter', 'latex')
ylabel('$F_z$', 'interpreter', 'latex')
