function [t_file t_actual] = cgns_read_flow_time(casename)

%% Reads the flow times of the cgns files given in the output
% file, as specified by flow-**.****.cgns etc (* is not
% indicative of sigfigs in this case, there can be more or less as
% specified in the record.config file

% t_file is the time written in the cgns file name, i.e. flow-12.31415.cgns
% t_actual is the actual time written in the cgns file, found using h5read

%% Get the path to the output directory, and find the cgns files
path = [casename '/output'];
od = cd(path);
contents = dir('*cgns');
names = {contents(:).name};
flownames = names(strncmp(names, 'flow', 4));

%% Remove strings from flow file names and sort by number
numstr = strrep(flownames, 'flow-', '');
numstr = strrep(numstr, '.cgns', '');
num = str2double(strrep(numstr, '.', ''));
[dummy, index] = sort(num);

%% Reconstruct the time given in the file names in time-order
t_file = numstr(index);
flownames = flownames(index);

%% Read the cgns files and get the actual time
t_actual = zeros(1, length(t_file));
t_loc = '/Base/Zone0/Etc/Time/ data';

for i = 1:length(t_file)
  t_path = [path '/' flownames{i}];
  t_actual(i) = h5read(t_path, t_loc);
end

cd(od);
