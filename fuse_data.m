function test_series = fuse_data( velmat, locationmat, nor )
%% function fuse_data
% function test_series = fuse_data( velmat, locationmat, ldmmat )
% 
% DESCRIPTION
% The function takes the input from LIDAR, GPS and LDM and fusion it into
% one matrix. Double time entries will be merged. The output matrix is
% sorted by time in ascending order. THe number of blade rotations are not
% mandatory. The resulting struct will be without column 11. 
%
% INPUT
%  - vel_mat: numerical array of presicion double. The following data is
%   stored:
%    column 1: range in m. Radial distance from LIDAR
%    column 2: Radial velocity in m/s
%    column 3: Doppler intensity
%    column 4: Time in seconds, starting at 12AM of the measurement day
%    column 5: azimuthal angle in degrees
%    column 6: elevation angle in degrees
%    column 7: pitch in ??
%    column 8: roll in ??
% - locationmat: numerical array of precision double
%    Column 1: Time in seconds, starting at 12AM of the day.
%    Column 2: Northing in m of LIDAR
%    Column 3: Easting in m of LIDAR
%  - nor: n x 2 numerical array.
%    Column 1: Array of the averaged times and
%    Column 2: The computed number of revolutions (nor) in seconds.
% 
% OUTPUT
% - test_series: 
%    Column 1: Time in seconds, starting at 12AM of the measurement day
%    column 2: range in m. Radial distance from LIDAR
%    column 3: Radial velocity in m/s
%    column 4: Doppler intensity
%    column 5: azimuthal angle in degrees
%    column 6: elevation angle in degrees
%    column 7: pitch in ??
%    column 8: roll in ??
%    Column 9: Northing in m of LIDAR
%    Column 10: Easting in m of LIDAR
%    Column 11: The computed number of revolutions (nor) in seconds. (if
%    present
%
% Code by: Markus Schmidt
%
% $Revision: 1.2$ $Date: 2013/04/25 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Input Check
if nargin == 3
    isnor = true;
elseif nargin == 2
    isnor = false;
else
    error('Wrong number of input arguments.')
end

% Compute length of new time array. Assuming no identical times so that
% every entry from input parametres requires a slot.
l_vel = length(velmat(:,4));
l_loc = length(locationmat(:,1));

if isnor
    l_nor = length(nor(:,1));
    time_length = l_vel + l_loc + l_nor;
else
    time_length = l_vel + l_loc;
end

% Creation of struct array for saving of data, according to the values of
% input parametres.
test_series.time    = NaN(time_length,1);
test_series.range   = NaN(time_length,1);
test_series.vel_r   = NaN(time_length,1);
test_series.doppler = NaN(time_length,1);
test_series.az      = NaN(time_length,1);
test_series.el      = NaN(time_length,1);
test_series.pitch   = NaN(time_length,1);
test_series.roll    = NaN(time_length,1);
test_series.lat     = NaN(time_length,1);
test_series.lon     = NaN(time_length,1);
if isnor
test_series.nor     = NaN(time_length,1);
end

% Write data into struct
% From velmat
test_series.time(1:l_vel)    = velmat(:,4);
test_series.range(1:l_vel)   = velmat(:,1);
test_series.vel_r(1:l_vel)   = velmat(:,2);
test_series.doppler(1:l_vel) = velmat(:,3);
test_series.az(1:l_vel)      = velmat(:,5);
test_series.el(1:l_vel)      = velmat(:,6);
test_series.pitch(1:l_vel)   = velmat(:,7);
test_series.roll(1:l_vel)    = velmat(:,8);

% From locationmat
test_series.time(l_vel+1 : l_vel + l_loc)	 = locationmat(:,1);
test_series.lat(l_vel+1 : l_vel + l_loc)	 = locationmat(:,2);
test_series.lon(l_vel+1 : l_vel + l_loc)	 = locationmat(:,3);

% From nor
if isnor
test_series.time(l_vel+l_loc+1 : end)	 = nor(:,1);
test_series.nor(l_vel+l_loc+1 : end)	 = nor(:,2);
end

% Convert struct tmatrix in order to compute and sort it
test_series = struct2cell(test_series);
test_series = cell2mat(test_series');

% Deleting double entries
test_series = unique(test_series,'rows');

end
