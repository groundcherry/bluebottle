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

function [t_flow_str t_flow] = cgns_read_flow_time(casename);

%% Reads the flow times of the cgns files given in the output
% file, as specified by part-**.****.cgns etc (* is not
% indicative of sigfigs in this case, there can be more or less as
% specified in the record.config file

% t_flow_str - time written in cgns filename, i.e. flow-3.14159.cgns
% t_flow - time written in cgns file, using h5read

path = [casename '/output'];
od = cd(path);
contents = dir;
% Delete the . and .. directories from the contents
contents(1) = [];
contents(1) = [];
% Sort the contents by date added and remove unnecessary info
S = [contents(:).datenum];
[S,S] = sort(S);
contents = {contents(S).name};

% initialize time struct
t_flow_str = struct([]);

% Initialize k which increments to the next spot in the struct
j = 1;

% Look through directory contents and place flow times in the
% correct structure
for i = 1:length(contents)
    if strncmp('flow', contents{i},4) == 1
        time = sscanf(contents{i}, 'flow-%s');
        time = strrep(time, '.cgns', '');
        t_flow_str(j).time = time;
        tsol = '/Base/Zone0/Etc/Time/ data';
        t_path = [path '/' contents{i}];
        t_flow(j) = h5read(t_path, tsol);
        j = j+1;
    end
end

cd(od);
