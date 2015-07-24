%% pull_part_data.m
% Purpose: pulls part data from a given set of cgns files and saves it as a 
%   .mat file
%
% pull_part_data(dir, ts, te, vargin)
% INPUTS:
%         dir -- the simulation directory you wish to work with
%         ts  -- the starting time to pull
%         te  -- the ending time to pull
%         varargin
%             -- 'append': append data to an existing .mat file
function pull_part_data(dir, ts, te, options)

fprintf('Initializing... \n');

% go through options
if nargin == 4
  switch options
    case 'append'
      % append data to pre-existing file
      fprintf('\t''append'' option enabled\n');
      load part_data.mat
      ts = time(end);
      if (te <= ts)
        error('te <= ts')
      end
    otherwise
      fprintf('Unrecognized option. Current options are:\n');
      fprintf('\t append');
      error('Correct function inputs');
  end
end

% Read part time
[tstr tnum] = cgns_read_part_time(dir);
nInd = 1:length(tnum);

% 'cut' the time array b/t ts and te values, inclusive
ind = find(tnum < ts | tnum > te);
nInd(ind) = [];
tnum(ind) = [];
ts = nInd(1);
te = nInd(end);

% number of time steps
nt = length(nInd);

% number of particles
[temp, ~, ~] = cgns_read_part_position(dir, 0);
np = size(temp,1);

if nargin == 4
  switch options
    case 'append'
      % extend old arrays
      Xp = [Xp, zeros(np, nt)];
      Yp = [Yp, zeros(np, nt)];
      Zp = [Zp, zeros(np, nt)];
      Up = [Up, zeros(np, nt)];
      Vp = [Vp, zeros(np, nt)];
      Wp = [Wp, zeros(np, nt)];
      FX = [FX, zeros(np, nt)];
      FY = [FY, zeros(np, nt)];
      FZ = [FZ, zeros(np, nt)];
      FXi = [FXi, zeros(np, nt)];
      FYi = [FYi, zeros(np, nt)];
      FZi = [FZi, zeros(np, nt)];
      FXh = [FXh, zeros(np, nt)];
      FYh = [FYh, zeros(np, nt)];
      FZh = [FZh, zeros(np, nt)];
      time = [time, tnum];
  end
else
  % create new arrays
  Xp = zeros(np, nt);
  Yp = zeros(np, nt);
  Zp = zeros(np, nt);
  Up = zeros(np, nt);
  Vp = zeros(np, nt);
  Wp = zeros(np, nt);
  FX = zeros(np, nt);
  FY = zeros(np, nt);
  FZ = zeros(np, nt);
  FXi = zeros(np, nt);
  FYi = zeros(np, nt);
  FZi = zeros(np, nt);
  FXh = zeros(np, nt);
  FYh = zeros(np, nt);
  FZh = zeros(np, nt);
  time = tnum;
end

fprintf('Reading data... ');
% read part variables
nmsg = 0;
count = 0;
for i = ts:te
  count = count + 1;
  [Xp(:,i), Yp(:,i), Zp(:,i)] = cgns_read_part_position(dir, tstr{nInd(count)});
  [Up(:,i), Vp(:,i), Wp(:,i)] = cgns_read_part_vel(dir, tstr{nInd(count)});
  [FX(:,i), FY(:,i), FZ(:,i)] = cgns_read_part_force_total(dir, ...
                                    tstr{nInd(count)});
  [FXi(:,i), FYi(:,i), FZi(:,i)] = cgns_read_part_force_interaction(dir, ...
                                    tstr{nInd(count)});
  [FXh(:,i), FYh(:,i), FZh(:,i)] = cgns_read_part_force_hydro(dir, ...
                                    tstr{nInd(count)});

  msg = sprintf('%d of %d', count, nt);
  fprintf(repmat('\b', 1, nmsg));
  fprintf(msg);
  nmsg = numel(msg);
end

save('part_data.mat', 'time',...
     'Xp', 'Yp', 'Zp', ...
     'Up', 'Vp', 'Wp', ...
     'FX', 'FY', 'FZ', ...
     'FXi', 'FYi', 'FZi', ...
     'FXh', 'FYh', 'FZh');

fprintf('... Done!\n');     
