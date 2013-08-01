function blade_props = blade_nor( ldm_data )
%% function blade_nor
% blade_props = blade_nor( ldm_data )
% 
% DESCRIPTION
% The function calculates the number of revolutions. Furthermore, it
% identifies outliers by deviation and replaces them with a linear
% interpolation of its left and right neighbour.
%
% INPUT
% - ldm_data: .txt file with information of blade passage from LDM. The
% time intervalls are given in miliseconds
% 
% OUTPUT
%  - blade_props: n x 2 numerical array.
%    Column 1: Array of the averaged times and
%    Column 2: The computed number of revolutions (nor) in rad/seconds.
%
% Code by: Markus Schmidt
%
% $Revision: 2.0$ $Date: 2013/04/25 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Error checks
if nargin~=1
    error('Wrong number of input arguments.')
end

% Import of file:
ldmmat = dlmread(ldm_data);

% Filter data
diff_ldm = diff(ldmmat);
diff_true = (diff_ldm <= 1.5*median(diff_ldm)) &...
    (diff_ldm >= 0.5*median(diff_ldm));

% Compute number of revolutions
nor = 2*pi/3./diff_ldm;

% Compute the averaged times
time_av = diff_ldm + ldmmat(1:end-1,1);

blade_props = [time_av, nor];
blade_props = blade_props(diff_true,:);

end
