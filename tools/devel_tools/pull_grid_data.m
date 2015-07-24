%% pull_grid_data.m
% Purpose: pulls grid data from a given set of cgns files and saves it as a 
%   .mat file
%
% pull_grid_data(dir, ts, te, vargin)
% INPUTS:
%         dir -- the simulation directory you wish to work with
function pull_grid_data(dir)

fprintf('Initializing... \n');

% number of particles
[temp, ~, ~] = cgns_read_part_position(dir, 0);
np = size(temp,1);

fprintf('Reading data... ');
% read grid variables
[x, y, z] = cgns_read_grid(pwd);
dom.xs = min(x(:,1,1));
dom.xe = max(x(:,1,1));
dom.xl = dom.xe - dom.xs;
dom.dx = mean(diff(x(:,1,1)));
dom.ys = min(y(1,:,1));
dom.ye = max(y(1,:,1));
dom.yl = dom.ye - dom.ys;
dom.dy = mean(diff(y(1,:,1)));
dom.zs = min(z(1,1,:));
dom.ze = max(z(1,1,:));
dom.zl = dom.ze - dom.zs;
dom.dz = mean(diff(z(1,1,:)));
dom.N = np;

save('grid_data.mat', 'dom');

fprintf('... Done!\n');     
