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

function coeffs = cgns_read_part_coeffs(casename, time)
% CGNS_READ_PART_RADIUS Read the particle Lamb's coefficients from a
%   BLUEBOTTLE-generated CGNS file.
%
%   coeffs = CGNS_READ_PART_RADIUS (CASENAME, TIME) reads the particle Lamb's
%   coefficients from the simulation CASENAME at time TIME. The result is a
%   struct containing the following members (n and m are replaced with their
%   appropriate non-negative integer values less than or equal to the truncation
%   order of the Lamb's expansion used):
%     * pnm_re
%     * pnm_im
%     * phinm_re
%     * phinm_im
%     * chinm_re
%     * chinm_im
%   Each coefficient member contains N values, where N is the number of
%   particles.
%
%   Example:
%     cgns_read_part_coeffs('simulation', 3.14159) will read the
%     appropriate output file located in 'simulation/output'

% determine input type of 'time'

if isa(time, 'double') == 1
    % find the directory contents
    path = [casename '/output'];
    od = cd(path);
    contents = dir;
    % Get names of files, remove first two (..,.), find flows/parts
    names = {contents.name}';
    names(1:2) = [];
    check = find(strncmp(names, 'part', 4));
    % Get sigfigs
    sigfigs = names{check(1)};
    sigfigs = length(sigfigs(6:end-5)) - 2;
    t_format = sprintf('%%1.%df', sigfigs);
    tt = sprintf(t_format, time);
    cd(od);
elseif isa(time, 'char') == 1
    tt = time;
end

path = [casename '/output/part-' tt '.cgns'];

osol = '/Base/Zone0/Solution/LambOrder/ data';

o = h5read(path, osol);
[n, m] = size(o);

switch o(1)
  case 0
      p00_r_sol =   '/Base/Zone0/Solution/p00_re/ data';
      p00_i_sol =   '/Base/Zone0/Solution/p00_im/ data';
    phi00_r_sol = '/Base/Zone0/Solution/phi00_re/ data';
    phi00_i_sol = '/Base/Zone0/Solution/phi00_im/ data';
    chi00_r_sol = '/Base/Zone0/Solution/chi00_re/ data';
    chi00_i_sol = '/Base/Zone0/Solution/chi00_im/ data';

      p00_r = h5read(path,   p00_r_sol);
      p00_i = h5read(path,   p00_i_sol);
    phi00_r = h5read(path, phi00_r_sol);
    phi00_i = h5read(path, phi00_i_sol);
    chi00_r = h5read(path, chi00_r_sol);
    chi00_i = h5read(path, chi00_i_sol);

    coeffs = struct(  'p00_re',   p00_r,...
                      'p00_im',   p00_i,...
                    'phi00_re', phi00_r,...
                    'phi00_im', phi00_i,...
                    'chi00_re', chi00_r,...
                    'chi00_im', chi00_i);
  case 1
      p00_r_sol =   '/Base/Zone0/Solution/p00_re/ data';
      p00_i_sol =   '/Base/Zone0/Solution/p00_im/ data';
    phi00_r_sol = '/Base/Zone0/Solution/phi00_re/ data';
    phi00_i_sol = '/Base/Zone0/Solution/phi00_im/ data';
    chi00_r_sol = '/Base/Zone0/Solution/chi00_re/ data';
    chi00_i_sol = '/Base/Zone0/Solution/chi00_im/ data';
      p10_r_sol =   '/Base/Zone0/Solution/p10_re/ data';
      p10_i_sol =   '/Base/Zone0/Solution/p10_im/ data';
    phi10_r_sol = '/Base/Zone0/Solution/phi10_re/ data';
    phi10_i_sol = '/Base/Zone0/Solution/phi10_im/ data';
    chi10_r_sol = '/Base/Zone0/Solution/chi10_re/ data';
    chi10_i_sol = '/Base/Zone0/Solution/chi10_im/ data';
      p11_r_sol =   '/Base/Zone0/Solution/p11_re/ data';
      p11_i_sol =   '/Base/Zone0/Solution/p11_im/ data';
    phi11_r_sol = '/Base/Zone0/Solution/phi11_re/ data';
    phi11_i_sol = '/Base/Zone0/Solution/phi11_im/ data';
    chi11_r_sol = '/Base/Zone0/Solution/chi11_re/ data';
    chi11_i_sol = '/Base/Zone0/Solution/chi11_im/ data';

      p00_r = h5read(path,   p00_r_sol);
      p00_i = h5read(path,   p00_i_sol);
    phi00_r = h5read(path, phi00_r_sol);
    phi00_i = h5read(path, phi00_i_sol);
    chi00_r = h5read(path, chi00_r_sol);
    chi00_i = h5read(path, chi00_i_sol);
      p10_r = h5read(path,   p10_r_sol);
      p10_i = h5read(path,   p10_i_sol);
    phi10_r = h5read(path, phi10_r_sol);
    phi10_i = h5read(path, phi10_i_sol);
    chi10_r = h5read(path, chi10_r_sol);
    chi10_i = h5read(path, chi10_i_sol);
      p11_r = h5read(path,   p11_r_sol);
      p11_i = h5read(path,   p11_i_sol);
    phi11_r = h5read(path, phi11_r_sol);
    phi11_i = h5read(path, phi11_i_sol);
    chi11_r = h5read(path, chi11_r_sol);
    chi11_i = h5read(path, chi11_i_sol);

    coeffs = struct(  'p00_re',   p00_r,...
                      'p00_im',   p00_i,...
                    'phi00_re', phi00_r,...
                    'phi00_im', phi00_i,...
                    'chi00_re', chi00_r,...
                    'chi00_im', chi00_i,...
                      'p10_re',   p10_r,...
                      'p10_im',   p10_i,...
                    'phi10_re', phi10_r,...
                    'phi10_im', phi10_i,...
                    'chi10_re', chi10_r,...
                    'chi10_im', chi10_i,...
                      'p11_re',   p11_r,...
                      'p11_im',   p11_i,...
                    'phi11_re', phi11_r,...
                    'phi11_im', phi11_i,...
                    'chi11_re', chi11_r,...
                    'chi11_im', chi11_i);
  case 2
      p00_r_sol =   '/Base/Zone0/Solution/p00_re/ data';
      p00_i_sol =   '/Base/Zone0/Solution/p00_im/ data';
    phi00_r_sol = '/Base/Zone0/Solution/phi00_re/ data';
    phi00_i_sol = '/Base/Zone0/Solution/phi00_im/ data';
    chi00_r_sol = '/Base/Zone0/Solution/chi00_re/ data';
    chi00_i_sol = '/Base/Zone0/Solution/chi00_im/ data';
      p10_r_sol =   '/Base/Zone0/Solution/p10_re/ data';
      p10_i_sol =   '/Base/Zone0/Solution/p10_im/ data';
    phi10_r_sol = '/Base/Zone0/Solution/phi10_re/ data';
    phi10_i_sol = '/Base/Zone0/Solution/phi10_im/ data';
    chi10_r_sol = '/Base/Zone0/Solution/chi10_re/ data';
    chi10_i_sol = '/Base/Zone0/Solution/chi10_im/ data';
      p11_r_sol =   '/Base/Zone0/Solution/p11_re/ data';
      p11_i_sol =   '/Base/Zone0/Solution/p11_im/ data';
    phi11_r_sol = '/Base/Zone0/Solution/phi11_re/ data';
    phi11_i_sol = '/Base/Zone0/Solution/phi11_im/ data';
    chi11_r_sol = '/Base/Zone0/Solution/chi11_re/ data';
    chi11_i_sol = '/Base/Zone0/Solution/chi11_im/ data';
      p20_r_sol =   '/Base/Zone0/Solution/p20_re/ data';
      p20_i_sol =   '/Base/Zone0/Solution/p20_im/ data';
    phi20_r_sol = '/Base/Zone0/Solution/phi20_re/ data';
    phi20_i_sol = '/Base/Zone0/Solution/phi20_im/ data';
    chi20_r_sol = '/Base/Zone0/Solution/chi20_re/ data';
    chi20_i_sol = '/Base/Zone0/Solution/chi20_im/ data';
      p21_r_sol =   '/Base/Zone0/Solution/p21_re/ data';
      p21_i_sol =   '/Base/Zone0/Solution/p21_im/ data';
    phi21_r_sol = '/Base/Zone0/Solution/phi21_re/ data';
    phi21_i_sol = '/Base/Zone0/Solution/phi21_im/ data';
    chi21_r_sol = '/Base/Zone0/Solution/chi21_re/ data';
    chi21_i_sol = '/Base/Zone0/Solution/chi21_im/ data';
      p22_r_sol =   '/Base/Zone0/Solution/p22_re/ data';
      p22_i_sol =   '/Base/Zone0/Solution/p22_im/ data';
    phi22_r_sol = '/Base/Zone0/Solution/phi22_re/ data';
    phi22_i_sol = '/Base/Zone0/Solution/phi22_im/ data';
    chi22_r_sol = '/Base/Zone0/Solution/chi22_re/ data';
    chi22_i_sol = '/Base/Zone0/Solution/chi22_im/ data';

      p00_r = h5read(path,   p00_r_sol);
      p00_i = h5read(path,   p00_i_sol);
    phi00_r = h5read(path, phi00_r_sol);
    phi00_i = h5read(path, phi00_i_sol);
    chi00_r = h5read(path, chi00_r_sol);
    chi00_i = h5read(path, chi00_i_sol);
      p10_r = h5read(path,   p10_r_sol);
      p10_i = h5read(path,   p10_i_sol);
    phi10_r = h5read(path, phi10_r_sol);
    phi10_i = h5read(path, phi10_i_sol);
    chi10_r = h5read(path, chi10_r_sol);
    chi10_i = h5read(path, chi10_i_sol);
      p11_r = h5read(path,   p11_r_sol);
      p11_i = h5read(path,   p11_i_sol);
    phi11_r = h5read(path, phi11_r_sol);
    phi11_i = h5read(path, phi11_i_sol);
    chi11_r = h5read(path, chi11_r_sol);
    chi11_i = h5read(path, chi11_i_sol);
      p20_r = h5read(path,   p20_r_sol);
      p20_i = h5read(path,   p20_i_sol);
    phi20_r = h5read(path, phi20_r_sol);
    phi20_i = h5read(path, phi20_i_sol);
    chi20_r = h5read(path, chi20_r_sol);
    chi20_i = h5read(path, chi20_i_sol);
      p21_r = h5read(path,   p21_r_sol);
      p21_i = h5read(path,   p21_i_sol);
    phi21_r = h5read(path, phi21_r_sol);
    phi21_i = h5read(path, phi21_i_sol);
    chi21_r = h5read(path, chi21_r_sol);
    chi21_i = h5read(path, chi21_i_sol);
      p22_r = h5read(path,   p22_r_sol);
      p22_i = h5read(path,   p22_i_sol);
    phi22_r = h5read(path, phi22_r_sol);
    phi22_i = h5read(path, phi22_i_sol);
    chi22_r = h5read(path, chi22_r_sol);
    chi22_i = h5read(path, chi22_i_sol);

    coeffs = struct(  'p00_re',   p00_r,...
                      'p00_im',   p00_i,...
                    'phi00_re', phi00_r,...
                    'phi00_im', phi00_i,...
                    'chi00_re', chi00_r,...
                    'chi00_im', chi00_i,...
                      'p10_re',   p10_r,...
                      'p10_im',   p10_i,...
                    'phi10_re', phi10_r,...
                    'phi10_im', phi10_i,...
                    'chi10_re', chi10_r,...
                    'chi10_im', chi10_i,...
                      'p11_re',   p11_r,...
                      'p11_im',   p11_i,...
                    'phi11_re', phi11_r,...
                    'phi11_im', phi11_i,...
                    'chi11_re', chi11_r,...
                    'chi11_im', chi11_i,...
                      'p20_re',   p20_r,...
                      'p20_im',   p20_i,...
                    'phi20_re', phi20_r,...
                    'phi20_im', phi20_i,...
                    'chi20_re', chi20_r,...
                    'chi20_im', chi20_i,...
                      'p21_re',   p21_r,...
                      'p21_im',   p21_i,...
                    'phi21_re', phi21_r,...
                    'phi21_im', phi21_i,...
                    'chi21_re', chi21_r,...
                    'chi21_im', chi21_i,...
                      'p22_re',   p22_r,...
                      'p22_im',   p22_i,...
                    'phi22_re', phi22_r,...
                    'phi22_im', phi22_i,...
                    'chi22_re', chi22_r,...
                    'chi22_im', chi22_i);
  case 3
      p00_r_sol =   '/Base/Zone0/Solution/p00_re/ data';
      p00_i_sol =   '/Base/Zone0/Solution/p00_im/ data';
    phi00_r_sol = '/Base/Zone0/Solution/phi00_re/ data';
    phi00_i_sol = '/Base/Zone0/Solution/phi00_im/ data';
    chi00_r_sol = '/Base/Zone0/Solution/chi00_re/ data';
    chi00_i_sol = '/Base/Zone0/Solution/chi00_im/ data';
      p10_r_sol =   '/Base/Zone0/Solution/p10_re/ data';
      p10_i_sol =   '/Base/Zone0/Solution/p10_im/ data';
    phi10_r_sol = '/Base/Zone0/Solution/phi10_re/ data';
    phi10_i_sol = '/Base/Zone0/Solution/phi10_im/ data';
    chi10_r_sol = '/Base/Zone0/Solution/chi10_re/ data';
    chi10_i_sol = '/Base/Zone0/Solution/chi10_im/ data';
      p11_r_sol =   '/Base/Zone0/Solution/p11_re/ data';
      p11_i_sol =   '/Base/Zone0/Solution/p11_im/ data';
    phi11_r_sol = '/Base/Zone0/Solution/phi11_re/ data';
    phi11_i_sol = '/Base/Zone0/Solution/phi11_im/ data';
    chi11_r_sol = '/Base/Zone0/Solution/chi11_re/ data';
    chi11_i_sol = '/Base/Zone0/Solution/chi11_im/ data';
      p20_r_sol =   '/Base/Zone0/Solution/p20_re/ data';
      p20_i_sol =   '/Base/Zone0/Solution/p20_im/ data';
    phi20_r_sol = '/Base/Zone0/Solution/phi20_re/ data';
    phi20_i_sol = '/Base/Zone0/Solution/phi20_im/ data';
    chi20_r_sol = '/Base/Zone0/Solution/chi20_re/ data';
    chi20_i_sol = '/Base/Zone0/Solution/chi20_im/ data';
      p21_r_sol =   '/Base/Zone0/Solution/p21_re/ data';
      p21_i_sol =   '/Base/Zone0/Solution/p21_im/ data';
    phi21_r_sol = '/Base/Zone0/Solution/phi21_re/ data';
    phi21_i_sol = '/Base/Zone0/Solution/phi21_im/ data';
    chi21_r_sol = '/Base/Zone0/Solution/chi21_re/ data';
    chi21_i_sol = '/Base/Zone0/Solution/chi21_im/ data';
      p22_r_sol =   '/Base/Zone0/Solution/p22_re/ data';
      p22_i_sol =   '/Base/Zone0/Solution/p22_im/ data';
    phi22_r_sol = '/Base/Zone0/Solution/phi22_re/ data';
    phi22_i_sol = '/Base/Zone0/Solution/phi22_im/ data';
    chi22_r_sol = '/Base/Zone0/Solution/chi22_re/ data';
    chi22_i_sol = '/Base/Zone0/Solution/chi22_im/ data';
      p30_r_sol =   '/Base/Zone0/Solution/p30_re/ data';
      p30_i_sol =   '/Base/Zone0/Solution/p30_im/ data';
    phi30_r_sol = '/Base/Zone0/Solution/phi30_re/ data';
    phi30_i_sol = '/Base/Zone0/Solution/phi30_im/ data';
    chi30_r_sol = '/Base/Zone0/Solution/chi30_re/ data';
    chi30_i_sol = '/Base/Zone0/Solution/chi30_im/ data';
      p31_r_sol =   '/Base/Zone0/Solution/p31_re/ data';
      p31_i_sol =   '/Base/Zone0/Solution/p31_im/ data';
    phi31_r_sol = '/Base/Zone0/Solution/phi31_re/ data';
    phi31_i_sol = '/Base/Zone0/Solution/phi31_im/ data';
    chi31_r_sol = '/Base/Zone0/Solution/chi31_re/ data';
    chi31_i_sol = '/Base/Zone0/Solution/chi31_im/ data';
      p32_r_sol =   '/Base/Zone0/Solution/p32_re/ data';
      p32_i_sol =   '/Base/Zone0/Solution/p32_im/ data';
    phi32_r_sol = '/Base/Zone0/Solution/phi32_re/ data';
    phi32_i_sol = '/Base/Zone0/Solution/phi32_im/ data';
    chi32_r_sol = '/Base/Zone0/Solution/chi32_re/ data';
    chi32_i_sol = '/Base/Zone0/Solution/chi32_im/ data';
      p33_r_sol =   '/Base/Zone0/Solution/p33_re/ data';
      p33_i_sol =   '/Base/Zone0/Solution/p33_im/ data';
    phi33_r_sol = '/Base/Zone0/Solution/phi33_re/ data';
    phi33_i_sol = '/Base/Zone0/Solution/phi33_im/ data';
    chi33_r_sol = '/Base/Zone0/Solution/chi33_re/ data';
    chi33_i_sol = '/Base/Zone0/Solution/chi33_im/ data';

      p00_r = h5read(path,   p00_r_sol);
      p00_i = h5read(path,   p00_i_sol);
    phi00_r = h5read(path, phi00_r_sol);
    phi00_i = h5read(path, phi00_i_sol);
    chi00_r = h5read(path, chi00_r_sol);
    chi00_i = h5read(path, chi00_i_sol);
      p10_r = h5read(path,   p10_r_sol);
      p10_i = h5read(path,   p10_i_sol);
    phi10_r = h5read(path, phi10_r_sol);
    phi10_i = h5read(path, phi10_i_sol);
    chi10_r = h5read(path, chi10_r_sol);
    chi10_i = h5read(path, chi10_i_sol);
      p11_r = h5read(path,   p11_r_sol);
      p11_i = h5read(path,   p11_i_sol);
    phi11_r = h5read(path, phi11_r_sol);
    phi11_i = h5read(path, phi11_i_sol);
    chi11_r = h5read(path, chi11_r_sol);
    chi11_i = h5read(path, chi11_i_sol);
      p20_r = h5read(path,   p20_r_sol);
      p20_i = h5read(path,   p20_i_sol);
    phi20_r = h5read(path, phi20_r_sol);
    phi20_i = h5read(path, phi20_i_sol);
    chi20_r = h5read(path, chi20_r_sol);
    chi20_i = h5read(path, chi20_i_sol);
      p21_r = h5read(path,   p21_r_sol);
      p21_i = h5read(path,   p21_i_sol);
    phi21_r = h5read(path, phi21_r_sol);
    phi21_i = h5read(path, phi21_i_sol);
    chi21_r = h5read(path, chi21_r_sol);
    chi21_i = h5read(path, chi21_i_sol);
      p22_r = h5read(path,   p22_r_sol);
      p22_i = h5read(path,   p22_i_sol);
    phi22_r = h5read(path, phi22_r_sol);
    phi22_i = h5read(path, phi22_i_sol);
    chi22_r = h5read(path, chi22_r_sol);
    chi22_i = h5read(path, chi22_i_sol);
      p30_r = h5read(path,   p30_r_sol);
      p30_i = h5read(path,   p30_i_sol);
    phi30_r = h5read(path, phi30_r_sol);
    phi30_i = h5read(path, phi30_i_sol);
    chi30_r = h5read(path, chi30_r_sol);
    chi30_i = h5read(path, chi30_i_sol);
      p31_r = h5read(path,   p31_r_sol);
      p31_i = h5read(path,   p31_i_sol);
    phi31_r = h5read(path, phi31_r_sol);
    phi31_i = h5read(path, phi31_i_sol);
    chi31_r = h5read(path, chi31_r_sol);
    chi31_i = h5read(path, chi31_i_sol);
      p32_r = h5read(path,   p32_r_sol);
      p32_i = h5read(path,   p32_i_sol);
    phi32_r = h5read(path, phi32_r_sol);
    phi32_i = h5read(path, phi32_i_sol);
    chi32_r = h5read(path, chi32_r_sol);
    chi32_i = h5read(path, chi32_i_sol);
      p33_r = h5read(path,   p33_r_sol);
      p33_i = h5read(path,   p33_i_sol);
    phi33_r = h5read(path, phi33_r_sol);
    phi33_i = h5read(path, phi33_i_sol);
    chi33_r = h5read(path, chi33_r_sol);
    chi33_i = h5read(path, chi33_i_sol);

    coeffs = struct(  'p00_re',   p00_r,...
                      'p00_im',   p00_i,...
                    'phi00_re', phi00_r,...
                    'phi00_im', phi00_i,...
                    'chi00_re', chi00_r,...
                    'chi00_im', chi00_i,...
                      'p10_re',   p10_r,...
                      'p10_im',   p10_i,...
                    'phi10_re', phi10_r,...
                    'phi10_im', phi10_i,...
                    'chi10_re', chi10_r,...
                    'chi10_im', chi10_i,...
                      'p11_re',   p11_r,...
                      'p11_im',   p11_i,...
                    'phi11_re', phi11_r,...
                    'phi11_im', phi11_i,...
                    'chi11_re', chi11_r,...
                    'chi11_im', chi11_i,...
                      'p20_re',   p20_r,...
                      'p20_im',   p20_i,...
                    'phi20_re', phi20_r,...
                    'phi20_im', phi20_i,...
                    'chi20_re', chi20_r,...
                    'chi20_im', chi20_i,...
                      'p21_re',   p21_r,...
                      'p21_im',   p21_i,...
                    'phi21_re', phi21_r,...
                    'phi21_im', phi21_i,...
                    'chi21_re', chi21_r,...
                    'chi21_im', chi21_i,...
                      'p22_re',   p22_r,...
                      'p22_im',   p22_i,...
                    'phi22_re', phi22_r,...
                    'phi22_im', phi22_i,...
                    'chi22_re', chi22_r,...
                    'chi22_im', chi22_i,...
                      'p30_re',   p30_r,...
                      'p30_im',   p30_i,...
                    'phi30_re', phi30_r,...
                    'phi30_im', phi30_i,...
                    'chi30_re', chi30_r,...
                    'chi30_im', chi30_i,...
                      'p31_re',   p31_r,...
                      'p31_im',   p31_i,...
                    'phi31_re', phi31_r,...
                    'phi31_im', phi31_i,...
                    'chi31_re', chi31_r,...
                    'chi31_im', chi31_i,...
                      'p32_re',   p32_r,...
                      'p32_im',   p32_i,...
                    'phi32_re', phi32_r,...
                    'phi32_im', phi32_i,...
                    'chi32_re', chi32_r,...
                    'chi32_im', chi32_i,...
                      'p33_re',   p33_r,...
                      'p33_im',   p33_i,...
                    'phi33_re', phi33_r,...
                    'phi33_im', phi33_i,...
                    'chi33_re', chi33_r,...
                    'chi33_im', chi33_i);
  case 4
      p00_r_sol =   '/Base/Zone0/Solution/p00_re/ data';
      p00_i_sol =   '/Base/Zone0/Solution/p00_im/ data';
    phi00_r_sol = '/Base/Zone0/Solution/phi00_re/ data';
    phi00_i_sol = '/Base/Zone0/Solution/phi00_im/ data';
    chi00_r_sol = '/Base/Zone0/Solution/chi00_re/ data';
    chi00_i_sol = '/Base/Zone0/Solution/chi00_im/ data';
      p10_r_sol =   '/Base/Zone0/Solution/p10_re/ data';
      p10_i_sol =   '/Base/Zone0/Solution/p10_im/ data';
    phi10_r_sol = '/Base/Zone0/Solution/phi10_re/ data';
    phi10_i_sol = '/Base/Zone0/Solution/phi10_im/ data';
    chi10_r_sol = '/Base/Zone0/Solution/chi10_re/ data';
    chi10_i_sol = '/Base/Zone0/Solution/chi10_im/ data';
      p11_r_sol =   '/Base/Zone0/Solution/p11_re/ data';
      p11_i_sol =   '/Base/Zone0/Solution/p11_im/ data';
    phi11_r_sol = '/Base/Zone0/Solution/phi11_re/ data';
    phi11_i_sol = '/Base/Zone0/Solution/phi11_im/ data';
    chi11_r_sol = '/Base/Zone0/Solution/chi11_re/ data';
    chi11_i_sol = '/Base/Zone0/Solution/chi11_im/ data';
      p20_r_sol =   '/Base/Zone0/Solution/p20_re/ data';
      p20_i_sol =   '/Base/Zone0/Solution/p20_im/ data';
    phi20_r_sol = '/Base/Zone0/Solution/phi20_re/ data';
    phi20_i_sol = '/Base/Zone0/Solution/phi20_im/ data';
    chi20_r_sol = '/Base/Zone0/Solution/chi20_re/ data';
    chi20_i_sol = '/Base/Zone0/Solution/chi20_im/ data';
      p21_r_sol =   '/Base/Zone0/Solution/p21_re/ data';
      p21_i_sol =   '/Base/Zone0/Solution/p21_im/ data';
    phi21_r_sol = '/Base/Zone0/Solution/phi21_re/ data';
    phi21_i_sol = '/Base/Zone0/Solution/phi21_im/ data';
    chi21_r_sol = '/Base/Zone0/Solution/chi21_re/ data';
    chi21_i_sol = '/Base/Zone0/Solution/chi21_im/ data';
      p22_r_sol =   '/Base/Zone0/Solution/p22_re/ data';
      p22_i_sol =   '/Base/Zone0/Solution/p22_im/ data';
    phi22_r_sol = '/Base/Zone0/Solution/phi22_re/ data';
    phi22_i_sol = '/Base/Zone0/Solution/phi22_im/ data';
    chi22_r_sol = '/Base/Zone0/Solution/chi22_re/ data';
    chi22_i_sol = '/Base/Zone0/Solution/chi22_im/ data';
      p30_r_sol =   '/Base/Zone0/Solution/p30_re/ data';
      p30_i_sol =   '/Base/Zone0/Solution/p30_im/ data';
    phi30_r_sol = '/Base/Zone0/Solution/phi30_re/ data';
    phi30_i_sol = '/Base/Zone0/Solution/phi30_im/ data';
    chi30_r_sol = '/Base/Zone0/Solution/chi30_re/ data';
    chi30_i_sol = '/Base/Zone0/Solution/chi30_im/ data';
      p31_r_sol =   '/Base/Zone0/Solution/p31_re/ data';
      p31_i_sol =   '/Base/Zone0/Solution/p31_im/ data';
    phi31_r_sol = '/Base/Zone0/Solution/phi31_re/ data';
    phi31_i_sol = '/Base/Zone0/Solution/phi31_im/ data';
    chi31_r_sol = '/Base/Zone0/Solution/chi31_re/ data';
    chi31_i_sol = '/Base/Zone0/Solution/chi31_im/ data';
      p32_r_sol =   '/Base/Zone0/Solution/p32_re/ data';
      p32_i_sol =   '/Base/Zone0/Solution/p32_im/ data';
    phi32_r_sol = '/Base/Zone0/Solution/phi32_re/ data';
    phi32_i_sol = '/Base/Zone0/Solution/phi32_im/ data';
    chi32_r_sol = '/Base/Zone0/Solution/chi32_re/ data';
    chi32_i_sol = '/Base/Zone0/Solution/chi32_im/ data';
      p33_r_sol =   '/Base/Zone0/Solution/p33_re/ data';
      p33_i_sol =   '/Base/Zone0/Solution/p33_im/ data';
    phi33_r_sol = '/Base/Zone0/Solution/phi33_re/ data';
    phi33_i_sol = '/Base/Zone0/Solution/phi33_im/ data';
    chi33_r_sol = '/Base/Zone0/Solution/chi33_re/ data';
    chi33_i_sol = '/Base/Zone0/Solution/chi33_im/ data';
      p40_r_sol =   '/Base/Zone0/Solution/p40_re/ data';
      p40_i_sol =   '/Base/Zone0/Solution/p40_im/ data';
    phi40_r_sol = '/Base/Zone0/Solution/phi40_re/ data';
    phi40_i_sol = '/Base/Zone0/Solution/phi40_im/ data';
    chi40_r_sol = '/Base/Zone0/Solution/chi40_re/ data';
    chi40_i_sol = '/Base/Zone0/Solution/chi40_im/ data';
      p41_r_sol =   '/Base/Zone0/Solution/p41_re/ data';
      p41_i_sol =   '/Base/Zone0/Solution/p41_im/ data';
    phi41_r_sol = '/Base/Zone0/Solution/phi41_re/ data';
    phi41_i_sol = '/Base/Zone0/Solution/phi41_im/ data';
    chi41_r_sol = '/Base/Zone0/Solution/chi41_re/ data';
    chi41_i_sol = '/Base/Zone0/Solution/chi41_im/ data';
      p42_r_sol =   '/Base/Zone0/Solution/p42_re/ data';
      p42_i_sol =   '/Base/Zone0/Solution/p42_im/ data';
    phi42_r_sol = '/Base/Zone0/Solution/phi42_re/ data';
    phi42_i_sol = '/Base/Zone0/Solution/phi42_im/ data';
    chi42_r_sol = '/Base/Zone0/Solution/chi42_re/ data';
    chi42_i_sol = '/Base/Zone0/Solution/chi42_im/ data';
      p43_r_sol =   '/Base/Zone0/Solution/p43_re/ data';
      p43_i_sol =   '/Base/Zone0/Solution/p43_im/ data';
    phi43_r_sol = '/Base/Zone0/Solution/phi43_re/ data';
    phi43_i_sol = '/Base/Zone0/Solution/phi43_im/ data';
    chi43_r_sol = '/Base/Zone0/Solution/chi43_re/ data';
    chi43_i_sol = '/Base/Zone0/Solution/chi43_im/ data';
      p44_r_sol =   '/Base/Zone0/Solution/p44_re/ data';
      p44_i_sol =   '/Base/Zone0/Solution/p44_im/ data';
    phi44_r_sol = '/Base/Zone0/Solution/phi44_re/ data';
    phi44_i_sol = '/Base/Zone0/Solution/phi44_im/ data';
    chi44_r_sol = '/Base/Zone0/Solution/chi44_re/ data';
    chi44_i_sol = '/Base/Zone0/Solution/chi44_im/ data';

      p00_r = h5read(path,   p00_r_sol);
      p00_i = h5read(path,   p00_i_sol);
    phi00_r = h5read(path, phi00_r_sol);
    phi00_i = h5read(path, phi00_i_sol);
    chi00_r = h5read(path, chi00_r_sol);
    chi00_i = h5read(path, chi00_i_sol);
      p10_r = h5read(path,   p10_r_sol);
      p10_i = h5read(path,   p10_i_sol);
    phi10_r = h5read(path, phi10_r_sol);
    phi10_i = h5read(path, phi10_i_sol);
    chi10_r = h5read(path, chi10_r_sol);
    chi10_i = h5read(path, chi10_i_sol);
      p11_r = h5read(path,   p11_r_sol);
      p11_i = h5read(path,   p11_i_sol);
    phi11_r = h5read(path, phi11_r_sol);
    phi11_i = h5read(path, phi11_i_sol);
    chi11_r = h5read(path, chi11_r_sol);
    chi11_i = h5read(path, chi11_i_sol);
      p20_r = h5read(path,   p20_r_sol);
      p20_i = h5read(path,   p20_i_sol);
    phi20_r = h5read(path, phi20_r_sol);
    phi20_i = h5read(path, phi20_i_sol);
    chi20_r = h5read(path, chi20_r_sol);
    chi20_i = h5read(path, chi20_i_sol);
      p21_r = h5read(path,   p21_r_sol);
      p21_i = h5read(path,   p21_i_sol);
    phi21_r = h5read(path, phi21_r_sol);
    phi21_i = h5read(path, phi21_i_sol);
    chi21_r = h5read(path, chi21_r_sol);
    chi21_i = h5read(path, chi21_i_sol);
      p22_r = h5read(path,   p22_r_sol);
      p22_i = h5read(path,   p22_i_sol);
    phi22_r = h5read(path, phi22_r_sol);
    phi22_i = h5read(path, phi22_i_sol);
    chi22_r = h5read(path, chi22_r_sol);
    chi22_i = h5read(path, chi22_i_sol);
      p30_r = h5read(path,   p30_r_sol);
      p30_i = h5read(path,   p30_i_sol);
    phi30_r = h5read(path, phi30_r_sol);
    phi30_i = h5read(path, phi30_i_sol);
    chi30_r = h5read(path, chi30_r_sol);
    chi30_i = h5read(path, chi30_i_sol);
      p31_r = h5read(path,   p31_r_sol);
      p31_i = h5read(path,   p31_i_sol);
    phi31_r = h5read(path, phi31_r_sol);
    phi31_i = h5read(path, phi31_i_sol);
    chi31_r = h5read(path, chi31_r_sol);
    chi31_i = h5read(path, chi31_i_sol);
      p32_r = h5read(path,   p32_r_sol);
      p32_i = h5read(path,   p32_i_sol);
    phi32_r = h5read(path, phi32_r_sol);
    phi32_i = h5read(path, phi32_i_sol);
    chi32_r = h5read(path, chi32_r_sol);
    chi32_i = h5read(path, chi32_i_sol);
      p33_r = h5read(path,   p33_r_sol);
      p33_i = h5read(path,   p33_i_sol);
    phi33_r = h5read(path, phi33_r_sol);
    phi33_i = h5read(path, phi33_i_sol);
    chi33_r = h5read(path, chi33_r_sol);
    chi33_i = h5read(path, chi33_i_sol);
      p40_r = h5read(path,   p40_r_sol);
      p40_i = h5read(path,   p40_i_sol);
    phi40_r = h5read(path, phi40_r_sol);
    phi40_i = h5read(path, phi40_i_sol);
    chi40_r = h5read(path, chi40_r_sol);
    chi40_i = h5read(path, chi40_i_sol);
      p41_r = h5read(path,   p41_r_sol);
      p41_i = h5read(path,   p41_i_sol);
    phi41_r = h5read(path, phi41_r_sol);
    phi41_i = h5read(path, phi41_i_sol);
    chi41_r = h5read(path, chi41_r_sol);
    chi41_i = h5read(path, chi41_i_sol);
      p42_r = h5read(path,   p42_r_sol);
      p42_i = h5read(path,   p42_i_sol);
    phi42_r = h5read(path, phi42_r_sol);
    phi42_i = h5read(path, phi42_i_sol);
    chi42_r = h5read(path, chi42_r_sol);
    chi42_i = h5read(path, chi42_i_sol);
      p43_r = h5read(path,   p43_r_sol);
      p43_i = h5read(path,   p43_i_sol);
    phi43_r = h5read(path, phi43_r_sol);
    phi43_i = h5read(path, phi43_i_sol);
    chi43_r = h5read(path, chi43_r_sol);
    chi43_i = h5read(path, chi43_i_sol);
      p44_r = h5read(path,   p44_r_sol);
      p44_i = h5read(path,   p44_i_sol);
    phi44_r = h5read(path, phi44_r_sol);
    phi44_i = h5read(path, phi44_i_sol);
    chi44_r = h5read(path, chi44_r_sol);
    chi44_i = h5read(path, chi44_i_sol);

    coeffs = struct(  'p00_re',   p00_r,...
                      'p00_im',   p00_i,...
                    'phi00_re', phi00_r,...
                    'phi00_im', phi00_i,...
                    'chi00_re', chi00_r,...
                    'chi00_im', chi00_i,...
                      'p10_re',   p10_r,...
                      'p10_im',   p10_i,...
                    'phi10_re', phi10_r,...
                    'phi10_im', phi10_i,...
                    'chi10_re', chi10_r,...
                    'chi10_im', chi10_i,...
                      'p11_re',   p11_r,...
                      'p11_im',   p11_i,...
                    'phi11_re', phi11_r,...
                    'phi11_im', phi11_i,...
                    'chi11_re', chi11_r,...
                    'chi11_im', chi11_i,...
                      'p20_re',   p20_r,...
                      'p20_im',   p20_i,...
                    'phi20_re', phi20_r,...
                    'phi20_im', phi20_i,...
                    'chi20_re', chi20_r,...
                    'chi20_im', chi20_i,...
                      'p21_re',   p21_r,...
                      'p21_im',   p21_i,...
                    'phi21_re', phi21_r,...
                    'phi21_im', phi21_i,...
                    'chi21_re', chi21_r,...
                    'chi21_im', chi21_i,...
                      'p22_re',   p22_r,...
                      'p22_im',   p22_i,...
                    'phi22_re', phi22_r,...
                    'phi22_im', phi22_i,...
                    'chi22_re', chi22_r,...
                    'chi22_im', chi22_i,...
                      'p30_re',   p30_r,...
                      'p30_im',   p30_i,...
                    'phi30_re', phi30_r,...
                    'phi30_im', phi30_i,...
                    'chi30_re', chi30_r,...
                    'chi30_im', chi30_i,...
                      'p31_re',   p31_r,...
                      'p31_im',   p31_i,...
                    'phi31_re', phi31_r,...
                    'phi31_im', phi31_i,...
                    'chi31_re', chi31_r,...
                    'chi31_im', chi31_i,...
                      'p32_re',   p32_r,...
                      'p32_im',   p32_i,...
                    'phi32_re', phi32_r,...
                    'phi32_im', phi32_i,...
                    'chi32_re', chi32_r,...
                    'chi32_im', chi32_i,...
                      'p33_re',   p33_r,...
                      'p33_im',   p33_i,...
                    'phi33_re', phi33_r,...
                    'phi33_im', phi33_i,...
                    'chi33_re', chi33_r,...
                    'chi33_im', chi33_i,...
                      'p40_re',   p40_r,...
                      'p40_im',   p40_i,...
                    'phi40_re', phi40_r,...
                    'phi40_im', phi40_i,...
                    'chi40_re', chi40_r,...
                    'chi40_im', chi40_i,...
                      'p41_re',   p41_r,...
                      'p41_im',   p41_i,...
                    'phi41_re', phi41_r,...
                    'phi41_im', phi41_i,...
                    'chi41_re', chi41_r,...
                    'chi41_im', chi41_i,...
                      'p42_re',   p42_r,...
                      'p42_im',   p42_i,...
                    'phi42_re', phi42_r,...
                    'phi42_im', phi42_i,...
                    'chi42_re', chi42_r,...
                    'chi42_im', chi42_i,...
                      'p43_re',   p43_r,...
                      'p43_im',   p43_i,...
                    'phi43_re', phi43_r,...
                    'phi43_im', phi43_i,...
                    'chi43_re', chi43_r,...
                    'chi43_im', chi43_i,...
                      'p44_re',   p44_r,...
                      'p44_im',   p44_i,...
                    'phi44_re', phi44_r,...
                    'phi44_im', phi44_i,...
                    'chi44_re', chi44_r,...
                    'chi44_im', chi44_i);
  case 5
      p00_r_sol =   '/Base/Zone0/Solution/p00_re/ data';
      p00_i_sol =   '/Base/Zone0/Solution/p00_im/ data';
    phi00_r_sol = '/Base/Zone0/Solution/phi00_re/ data';
    phi00_i_sol = '/Base/Zone0/Solution/phi00_im/ data';
    chi00_r_sol = '/Base/Zone0/Solution/chi00_re/ data';
    chi00_i_sol = '/Base/Zone0/Solution/chi00_im/ data';
      p10_r_sol =   '/Base/Zone0/Solution/p10_re/ data';
      p10_i_sol =   '/Base/Zone0/Solution/p10_im/ data';
    phi10_r_sol = '/Base/Zone0/Solution/phi10_re/ data';
    phi10_i_sol = '/Base/Zone0/Solution/phi10_im/ data';
    chi10_r_sol = '/Base/Zone0/Solution/chi10_re/ data';
    chi10_i_sol = '/Base/Zone0/Solution/chi10_im/ data';
      p11_r_sol =   '/Base/Zone0/Solution/p11_re/ data';
      p11_i_sol =   '/Base/Zone0/Solution/p11_im/ data';
    phi11_r_sol = '/Base/Zone0/Solution/phi11_re/ data';
    phi11_i_sol = '/Base/Zone0/Solution/phi11_im/ data';
    chi11_r_sol = '/Base/Zone0/Solution/chi11_re/ data';
    chi11_i_sol = '/Base/Zone0/Solution/chi11_im/ data';
      p20_r_sol =   '/Base/Zone0/Solution/p20_re/ data';
      p20_i_sol =   '/Base/Zone0/Solution/p20_im/ data';
    phi20_r_sol = '/Base/Zone0/Solution/phi20_re/ data';
    phi20_i_sol = '/Base/Zone0/Solution/phi20_im/ data';
    chi20_r_sol = '/Base/Zone0/Solution/chi20_re/ data';
    chi20_i_sol = '/Base/Zone0/Solution/chi20_im/ data';
      p21_r_sol =   '/Base/Zone0/Solution/p21_re/ data';
      p21_i_sol =   '/Base/Zone0/Solution/p21_im/ data';
    phi21_r_sol = '/Base/Zone0/Solution/phi21_re/ data';
    phi21_i_sol = '/Base/Zone0/Solution/phi21_im/ data';
    chi21_r_sol = '/Base/Zone0/Solution/chi21_re/ data';
    chi21_i_sol = '/Base/Zone0/Solution/chi21_im/ data';
      p22_r_sol =   '/Base/Zone0/Solution/p22_re/ data';
      p22_i_sol =   '/Base/Zone0/Solution/p22_im/ data';
    phi22_r_sol = '/Base/Zone0/Solution/phi22_re/ data';
    phi22_i_sol = '/Base/Zone0/Solution/phi22_im/ data';
    chi22_r_sol = '/Base/Zone0/Solution/chi22_re/ data';
    chi22_i_sol = '/Base/Zone0/Solution/chi22_im/ data';
      p30_r_sol =   '/Base/Zone0/Solution/p30_re/ data';
      p30_i_sol =   '/Base/Zone0/Solution/p30_im/ data';
    phi30_r_sol = '/Base/Zone0/Solution/phi30_re/ data';
    phi30_i_sol = '/Base/Zone0/Solution/phi30_im/ data';
    chi30_r_sol = '/Base/Zone0/Solution/chi30_re/ data';
    chi30_i_sol = '/Base/Zone0/Solution/chi30_im/ data';
      p31_r_sol =   '/Base/Zone0/Solution/p31_re/ data';
      p31_i_sol =   '/Base/Zone0/Solution/p31_im/ data';
    phi31_r_sol = '/Base/Zone0/Solution/phi31_re/ data';
    phi31_i_sol = '/Base/Zone0/Solution/phi31_im/ data';
    chi31_r_sol = '/Base/Zone0/Solution/chi31_re/ data';
    chi31_i_sol = '/Base/Zone0/Solution/chi31_im/ data';
      p32_r_sol =   '/Base/Zone0/Solution/p32_re/ data';
      p32_i_sol =   '/Base/Zone0/Solution/p32_im/ data';
    phi32_r_sol = '/Base/Zone0/Solution/phi32_re/ data';
    phi32_i_sol = '/Base/Zone0/Solution/phi32_im/ data';
    chi32_r_sol = '/Base/Zone0/Solution/chi32_re/ data';
    chi32_i_sol = '/Base/Zone0/Solution/chi32_im/ data';
      p33_r_sol =   '/Base/Zone0/Solution/p33_re/ data';
      p33_i_sol =   '/Base/Zone0/Solution/p33_im/ data';
    phi33_r_sol = '/Base/Zone0/Solution/phi33_re/ data';
    phi33_i_sol = '/Base/Zone0/Solution/phi33_im/ data';
    chi33_r_sol = '/Base/Zone0/Solution/chi33_re/ data';
    chi33_i_sol = '/Base/Zone0/Solution/chi33_im/ data';
      p40_r_sol =   '/Base/Zone0/Solution/p40_re/ data';
      p40_i_sol =   '/Base/Zone0/Solution/p40_im/ data';
    phi40_r_sol = '/Base/Zone0/Solution/phi40_re/ data';
    phi40_i_sol = '/Base/Zone0/Solution/phi40_im/ data';
    chi40_r_sol = '/Base/Zone0/Solution/chi40_re/ data';
    chi40_i_sol = '/Base/Zone0/Solution/chi40_im/ data';
      p41_r_sol =   '/Base/Zone0/Solution/p41_re/ data';
      p41_i_sol =   '/Base/Zone0/Solution/p41_im/ data';
    phi41_r_sol = '/Base/Zone0/Solution/phi41_re/ data';
    phi41_i_sol = '/Base/Zone0/Solution/phi41_im/ data';
    chi41_r_sol = '/Base/Zone0/Solution/chi41_re/ data';
    chi41_i_sol = '/Base/Zone0/Solution/chi41_im/ data';
      p42_r_sol =   '/Base/Zone0/Solution/p42_re/ data';
      p42_i_sol =   '/Base/Zone0/Solution/p42_im/ data';
    phi42_r_sol = '/Base/Zone0/Solution/phi42_re/ data';
    phi42_i_sol = '/Base/Zone0/Solution/phi42_im/ data';
    chi42_r_sol = '/Base/Zone0/Solution/chi42_re/ data';
    chi42_i_sol = '/Base/Zone0/Solution/chi42_im/ data';
      p43_r_sol =   '/Base/Zone0/Solution/p43_re/ data';
      p43_i_sol =   '/Base/Zone0/Solution/p43_im/ data';
    phi43_r_sol = '/Base/Zone0/Solution/phi43_re/ data';
    phi43_i_sol = '/Base/Zone0/Solution/phi43_im/ data';
    chi43_r_sol = '/Base/Zone0/Solution/chi43_re/ data';
    chi43_i_sol = '/Base/Zone0/Solution/chi43_im/ data';
      p44_r_sol =   '/Base/Zone0/Solution/p44_re/ data';
      p44_i_sol =   '/Base/Zone0/Solution/p44_im/ data';
    phi44_r_sol = '/Base/Zone0/Solution/phi44_re/ data';
    phi44_i_sol = '/Base/Zone0/Solution/phi44_im/ data';
    chi44_r_sol = '/Base/Zone0/Solution/chi44_re/ data';
    chi44_i_sol = '/Base/Zone0/Solution/chi44_im/ data';
      p50_r_sol =   '/Base/Zone0/Solution/p50_re/ data';
      p50_i_sol =   '/Base/Zone0/Solution/p50_im/ data';
    phi50_r_sol = '/Base/Zone0/Solution/phi50_re/ data';
    phi50_i_sol = '/Base/Zone0/Solution/phi50_im/ data';
    chi50_r_sol = '/Base/Zone0/Solution/chi50_re/ data';
    chi50_i_sol = '/Base/Zone0/Solution/chi50_im/ data';
      p51_r_sol =   '/Base/Zone0/Solution/p51_re/ data';
      p51_i_sol =   '/Base/Zone0/Solution/p51_im/ data';
    phi51_r_sol = '/Base/Zone0/Solution/phi51_re/ data';
    phi51_i_sol = '/Base/Zone0/Solution/phi51_im/ data';
    chi51_r_sol = '/Base/Zone0/Solution/chi51_re/ data';
    chi51_i_sol = '/Base/Zone0/Solution/chi51_im/ data';
      p52_r_sol =   '/Base/Zone0/Solution/p52_re/ data';
      p52_i_sol =   '/Base/Zone0/Solution/p52_im/ data';
    phi52_r_sol = '/Base/Zone0/Solution/phi52_re/ data';
    phi52_i_sol = '/Base/Zone0/Solution/phi52_im/ data';
    chi52_r_sol = '/Base/Zone0/Solution/chi52_re/ data';
    chi52_i_sol = '/Base/Zone0/Solution/chi52_im/ data';
      p53_r_sol =   '/Base/Zone0/Solution/p53_re/ data';
      p53_i_sol =   '/Base/Zone0/Solution/p53_im/ data';
    phi53_r_sol = '/Base/Zone0/Solution/phi53_re/ data';
    phi53_i_sol = '/Base/Zone0/Solution/phi53_im/ data';
    chi53_r_sol = '/Base/Zone0/Solution/chi53_re/ data';
    chi53_i_sol = '/Base/Zone0/Solution/chi53_im/ data';
      p54_r_sol =   '/Base/Zone0/Solution/p54_re/ data';
      p54_i_sol =   '/Base/Zone0/Solution/p54_im/ data';
    phi54_r_sol = '/Base/Zone0/Solution/phi54_re/ data';
    phi54_i_sol = '/Base/Zone0/Solution/phi54_im/ data';
    chi54_r_sol = '/Base/Zone0/Solution/chi54_re/ data';
    chi54_i_sol = '/Base/Zone0/Solution/chi54_im/ data';
      p55_r_sol =   '/Base/Zone0/Solution/p55_re/ data';
      p55_i_sol =   '/Base/Zone0/Solution/p55_im/ data';
    phi55_r_sol = '/Base/Zone0/Solution/phi55_re/ data';
    phi55_i_sol = '/Base/Zone0/Solution/phi55_im/ data';
    chi55_r_sol = '/Base/Zone0/Solution/chi55_re/ data';
    chi55_i_sol = '/Base/Zone0/Solution/chi55_im/ data';

      p00_r = h5read(path,   p00_r_sol);
      p00_i = h5read(path,   p00_i_sol);
    phi00_r = h5read(path, phi00_r_sol);
    phi00_i = h5read(path, phi00_i_sol);
    chi00_r = h5read(path, chi00_r_sol);
    chi00_i = h5read(path, chi00_i_sol);
      p10_r = h5read(path,   p10_r_sol);
      p10_i = h5read(path,   p10_i_sol);
    phi10_r = h5read(path, phi10_r_sol);
    phi10_i = h5read(path, phi10_i_sol);
    chi10_r = h5read(path, chi10_r_sol);
    chi10_i = h5read(path, chi10_i_sol);
      p11_r = h5read(path,   p11_r_sol);
      p11_i = h5read(path,   p11_i_sol);
    phi11_r = h5read(path, phi11_r_sol);
    phi11_i = h5read(path, phi11_i_sol);
    chi11_r = h5read(path, chi11_r_sol);
    chi11_i = h5read(path, chi11_i_sol);
      p20_r = h5read(path,   p20_r_sol);
      p20_i = h5read(path,   p20_i_sol);
    phi20_r = h5read(path, phi20_r_sol);
    phi20_i = h5read(path, phi20_i_sol);
    chi20_r = h5read(path, chi20_r_sol);
    chi20_i = h5read(path, chi20_i_sol);
      p21_r = h5read(path,   p21_r_sol);
      p21_i = h5read(path,   p21_i_sol);
    phi21_r = h5read(path, phi21_r_sol);
    phi21_i = h5read(path, phi21_i_sol);
    chi21_r = h5read(path, chi21_r_sol);
    chi21_i = h5read(path, chi21_i_sol);
      p22_r = h5read(path,   p22_r_sol);
      p22_i = h5read(path,   p22_i_sol);
    phi22_r = h5read(path, phi22_r_sol);
    phi22_i = h5read(path, phi22_i_sol);
    chi22_r = h5read(path, chi22_r_sol);
    chi22_i = h5read(path, chi22_i_sol);
      p30_r = h5read(path,   p30_r_sol);
      p30_i = h5read(path,   p30_i_sol);
    phi30_r = h5read(path, phi30_r_sol);
    phi30_i = h5read(path, phi30_i_sol);
    chi30_r = h5read(path, chi30_r_sol);
    chi30_i = h5read(path, chi30_i_sol);
      p31_r = h5read(path,   p31_r_sol);
      p31_i = h5read(path,   p31_i_sol);
    phi31_r = h5read(path, phi31_r_sol);
    phi31_i = h5read(path, phi31_i_sol);
    chi31_r = h5read(path, chi31_r_sol);
    chi31_i = h5read(path, chi31_i_sol);
      p32_r = h5read(path,   p32_r_sol);
      p32_i = h5read(path,   p32_i_sol);
    phi32_r = h5read(path, phi32_r_sol);
    phi32_i = h5read(path, phi32_i_sol);
    chi32_r = h5read(path, chi32_r_sol);
    chi32_i = h5read(path, chi32_i_sol);
      p33_r = h5read(path,   p33_r_sol);
      p33_i = h5read(path,   p33_i_sol);
    phi33_r = h5read(path, phi33_r_sol);
    phi33_i = h5read(path, phi33_i_sol);
    chi33_r = h5read(path, chi33_r_sol);
    chi33_i = h5read(path, chi33_i_sol);
      p40_r = h5read(path,   p40_r_sol);
      p40_i = h5read(path,   p40_i_sol);
    phi40_r = h5read(path, phi40_r_sol);
    phi40_i = h5read(path, phi40_i_sol);
    chi40_r = h5read(path, chi40_r_sol);
    chi40_i = h5read(path, chi40_i_sol);
      p41_r = h5read(path,   p41_r_sol);
      p41_i = h5read(path,   p41_i_sol);
    phi41_r = h5read(path, phi41_r_sol);
    phi41_i = h5read(path, phi41_i_sol);
    chi41_r = h5read(path, chi41_r_sol);
    chi41_i = h5read(path, chi41_i_sol);
      p42_r = h5read(path,   p42_r_sol);
      p42_i = h5read(path,   p42_i_sol);
    phi42_r = h5read(path, phi42_r_sol);
    phi42_i = h5read(path, phi42_i_sol);
    chi42_r = h5read(path, chi42_r_sol);
    chi42_i = h5read(path, chi42_i_sol);
      p43_r = h5read(path,   p43_r_sol);
      p43_i = h5read(path,   p43_i_sol);
    phi43_r = h5read(path, phi43_r_sol);
    phi43_i = h5read(path, phi43_i_sol);
    chi43_r = h5read(path, chi43_r_sol);
    chi43_i = h5read(path, chi43_i_sol);
      p44_r = h5read(path,   p44_r_sol);
      p44_i = h5read(path,   p44_i_sol);
    phi44_r = h5read(path, phi44_r_sol);
    phi44_i = h5read(path, phi44_i_sol);
    chi44_r = h5read(path, chi44_r_sol);
    chi44_i = h5read(path, chi44_i_sol);
      p50_r = h5read(path,   p50_r_sol);
      p50_i = h5read(path,   p50_i_sol);
    phi50_r = h5read(path, phi50_r_sol);
    phi50_i = h5read(path, phi50_i_sol);
    chi50_r = h5read(path, chi50_r_sol);
    chi50_i = h5read(path, chi50_i_sol);
      p51_r = h5read(path,   p51_r_sol);
      p51_i = h5read(path,   p51_i_sol);
    phi51_r = h5read(path, phi51_r_sol);
    phi51_i = h5read(path, phi51_i_sol);
    chi51_r = h5read(path, chi51_r_sol);
    chi51_i = h5read(path, chi51_i_sol);
      p52_r = h5read(path,   p52_r_sol);
      p52_i = h5read(path,   p52_i_sol);
    phi52_r = h5read(path, phi52_r_sol);
    phi52_i = h5read(path, phi52_i_sol);
    chi52_r = h5read(path, chi52_r_sol);
    chi52_i = h5read(path, chi52_i_sol);
      p53_r = h5read(path,   p53_r_sol);
      p53_i = h5read(path,   p53_i_sol);
    phi53_r = h5read(path, phi53_r_sol);
    phi53_i = h5read(path, phi53_i_sol);
    chi53_r = h5read(path, chi53_r_sol);
    chi53_i = h5read(path, chi53_i_sol);
      p54_r = h5read(path,   p54_r_sol);
      p54_i = h5read(path,   p54_i_sol);
    phi54_r = h5read(path, phi54_r_sol);
    phi54_i = h5read(path, phi54_i_sol);
    chi54_r = h5read(path, chi54_r_sol);
    chi54_i = h5read(path, chi54_i_sol);
      p55_r = h5read(path,   p55_r_sol);
      p55_i = h5read(path,   p55_i_sol);
    phi55_r = h5read(path, phi55_r_sol);
    phi55_i = h5read(path, phi55_i_sol);
    chi55_r = h5read(path, chi55_r_sol);
    chi55_i = h5read(path, chi55_i_sol);

    coeffs = struct(  'p00_re',   p00_r,...
                      'p00_im',   p00_i,...
                    'phi00_re', phi00_r,...
                    'phi00_im', phi00_i,...
                    'chi00_re', chi00_r,...
                    'chi00_im', chi00_i,...
                      'p10_re',   p10_r,...
                      'p10_im',   p10_i,...
                    'phi10_re', phi10_r,...
                    'phi10_im', phi10_i,...
                    'chi10_re', chi10_r,...
                    'chi10_im', chi10_i,...
                      'p11_re',   p11_r,...
                      'p11_im',   p11_i,...
                    'phi11_re', phi11_r,...
                    'phi11_im', phi11_i,...
                    'chi11_re', chi11_r,...
                    'chi11_im', chi11_i,...
                      'p20_re',   p20_r,...
                      'p20_im',   p20_i,...
                    'phi20_re', phi20_r,...
                    'phi20_im', phi20_i,...
                    'chi20_re', chi20_r,...
                    'chi20_im', chi20_i,...
                      'p21_re',   p21_r,...
                      'p21_im',   p21_i,...
                    'phi21_re', phi21_r,...
                    'phi21_im', phi21_i,...
                    'chi21_re', chi21_r,...
                    'chi21_im', chi21_i,...
                      'p22_re',   p22_r,...
                      'p22_im',   p22_i,...
                    'phi22_re', phi22_r,...
                    'phi22_im', phi22_i,...
                    'chi22_re', chi22_r,...
                    'chi22_im', chi22_i,...
                      'p30_re',   p30_r,...
                      'p30_im',   p30_i,...
                    'phi30_re', phi30_r,...
                    'phi30_im', phi30_i,...
                    'chi30_re', chi30_r,...
                    'chi30_im', chi30_i,...
                      'p31_re',   p31_r,...
                      'p31_im',   p31_i,...
                    'phi31_re', phi31_r,...
                    'phi31_im', phi31_i,...
                    'chi31_re', chi31_r,...
                    'chi31_im', chi31_i,...
                      'p32_re',   p32_r,...
                      'p32_im',   p32_i,...
                    'phi32_re', phi32_r,...
                    'phi32_im', phi32_i,...
                    'chi32_re', chi32_r,...
                    'chi32_im', chi32_i,...
                      'p33_re',   p33_r,...
                      'p33_im',   p33_i,...
                    'phi33_re', phi33_r,...
                    'phi33_im', phi33_i,...
                    'chi33_re', chi33_r,...
                    'chi33_im', chi33_i,...
                      'p40_re',   p40_r,...
                      'p40_im',   p40_i,...
                    'phi40_re', phi40_r,...
                    'phi40_im', phi40_i,...
                    'chi40_re', chi40_r,...
                    'chi40_im', chi40_i,...
                      'p41_re',   p41_r,...
                      'p41_im',   p41_i,...
                    'phi41_re', phi41_r,...
                    'phi41_im', phi41_i,...
                    'chi41_re', chi41_r,...
                    'chi41_im', chi41_i,...
                      'p42_re',   p42_r,...
                      'p42_im',   p42_i,...
                    'phi42_re', phi42_r,...
                    'phi42_im', phi42_i,...
                    'chi42_re', chi42_r,...
                    'chi42_im', chi42_i,...
                      'p43_re',   p43_r,...
                      'p43_im',   p43_i,...
                    'phi43_re', phi43_r,...
                    'phi43_im', phi43_i,...
                    'chi43_re', chi43_r,...
                    'chi43_im', chi43_i,...
                      'p44_re',   p44_r,...
                      'p44_im',   p44_i,...
                    'phi44_re', phi44_r,...
                    'phi44_im', phi44_i,...
                    'chi44_re', chi44_r,...
                    'chi44_im', chi44_i,...
                      'p50_re',   p50_r,...
                      'p50_im',   p50_i,...
                    'phi50_re', phi50_r,...
                    'phi50_im', phi50_i,...
                    'chi50_re', chi50_r,...
                    'chi50_im', chi50_i,...
                      'p51_re',   p51_r,...
                      'p51_im',   p51_i,...
                    'phi51_re', phi51_r,...
                    'phi51_im', phi51_i,...
                    'chi51_re', chi51_r,...
                    'chi51_im', chi51_i,...
                      'p52_re',   p52_r,...
                      'p52_im',   p52_i,...
                    'phi52_re', phi52_r,...
                    'phi52_im', phi52_i,...
                    'chi52_re', chi52_r,...
                    'chi52_im', chi52_i,...
                      'p53_re',   p53_r,...
                      'p53_im',   p53_i,...
                    'phi53_re', phi53_r,...
                    'phi53_im', phi53_i,...
                    'chi53_re', chi53_r,...
                    'chi53_im', chi53_i,...
                      'p54_re',   p54_r,...
                      'p54_im',   p54_i,...
                    'phi54_re', phi54_r,...
                    'phi54_im', phi54_i,...
                    'chi54_re', chi54_r,...
                    'chi54_im', chi54_i,...
                      'p55_re',   p55_r,...
                      'p55_im',   p55_i,...
                    'phi55_re', phi55_r,...
                    'phi55_im', phi55_i,...
                    'chi55_re', chi55_r,...
                    'chi55_im', chi55_i);
end
