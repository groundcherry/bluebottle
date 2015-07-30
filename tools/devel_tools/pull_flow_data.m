%% pull_flow_data.m
% Purpose: pulls flow data from a given set of cgns files and saves it as a 
%   .mat file
%
% pull_flow_data(dir, ts, te, vargin)
% INPUTS:
%         dir -- the simulation directory you wish to work with
%         ts  -- the starting time to pull
%         te  -- the ending time to pull
%         varargin
%             -- 'append': append data to an existing .mat file
function pull_flow_data(dir, ts, te, options)

fprintf('Initializing... \n');

% go through options
if nargin == 4
  switch vargin
    case 'append'
      % append data to pre-existing file
      fprintf('\t''append'' option enabled\n');
      load flow_data.mat
      ts = time(end);
      if (te <= ts)
        error('te <= ts')
      end
    otherwise
      fprintf('Unrecognized option. Current options are:\n');
      fprintf('\t append');
      error('Correct function inputs.')
  end
end

% Read flow time
[tstr tnum] = cgns_read_flow_time(dir);
nInd = 1:length(tnum);

% 'cut' the time array b/t ts and te values, inclusive
ind = find(tnum < ts | tnum > te);
nInd(ind) = [];
tnum(ind) = [];
ts = nInd(1);
te = nInd(end);

% number of time steps
nt = length(nInd);

% find dimensions of flow arrays
% -- everything interpolated to cell center, so one set of dims should be good
temp = cgns_read_flow_vel(pwd, tstr{nInd(1)});
[ni nj nk] = size(temp);

if nargin == 4
  switch options
    case 'append'
      % extend old arrays
      Uf = [Uf, zeros(ni, nj, nk, nt)];
      Vf = [Vf, zeros(ni, nj, nk, nt)];
      Wf = [Wf, zeros(ni, nj, nk, nt)];
      phase = [phase, zeros(ni, nj, nk, nt)];
      time = [time, tnum];
  end
else
  % create new arrays
  Uf = zeros(ni, nj, nk, nt);
  Vf = zeros(ni, nj, nk, nt);
  Wf = zeros(ni, nj, nk, nt);
  phase = zeros(ni, nj, nk, nt);
  time = tnum;
end



fprintf('Reading data... ');
% read part variables
nmsg = 0;
count = 0;
for i = ts:te
  count = count + 1;
  % read vel
  [Uf(:,:,:,i), Vf(:,:,:,i), Wf(:,:,:,i)] = cgns_read_flow_vel(dir, ...
                                                tstr{nInd(count)});
  % read phase
  phase(:,:,:,i) = cgns_read_flow_phase(dir, tstr{nInd(count)});

  msg = sprintf('%d of %d', count, nt);
  fprintf(repmat('\b', 1, nmsg));
  fprintf(msg);
  nmsg = numel(msg);
end

save('flow_data.mat', 'time', 'Uf', 'Vf', 'Wf', 'phase')

fprintf('... Done!\n');
