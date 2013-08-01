function [ test_series ] = UTM_interp(velmat, loc_props)
%% function UTM_interp
% [ test_series ] = UTM_interp(velmat, loc_props)
% 
% DESCRIPTION The function takes existing matrices velmat and loc_props. It
% interpolates the location to the existing measurements of lidar and
% deletes the exiting entries for the location. The interpolation is
% computed with 1D linear interpolation.
%
% INPUT
% - vel_mat: numerical array of presicion double. The following data is
% stored:
% column 1: range in m. Radial distance from LIDAR
% column 2: Radial velocity in m/s
% column 3: Doppler intensity
% column 4: Time in seconds, starting at 12AM of the measurement day
% column 5: azimuthal angle in degrees
% column 6: elevation angle in degrees
% column 7: pitch in ??
% column 8: roll in ??
% - loc_props: numerical array of precision double
%   Column 1: Time in seconds, starting at 12AM of the day.
%   Column 2: Northing in m of LIDAR
%   Column 3: Easting in m of LIDAR
%   Column 4: heading in rad of LIDAR (eastwards defines 0 rad)
%   Column 5: Velocity in m/s of LIDAR
%
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
%    column 9: Northing in m of LIDAR
%    column 10: Easting in m of LIDAR
%    column 11: heading in rad of LIDAR (eastwards defines 0 rad)
%    column 12: Velocity in m/s of LIDAR
%
% Code by: Markus Schmidt
%
% $Revision: 0.4$ $Date: 2013/05/14 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Input check
if nargin < 2 || nargin > 4
    error('Wrong number of input arguments')
end

% Merge data
test_series = NaN(size(velmat,1)+size(loc_props,1),12);

test_series(1:size(velmat,1),[2 3 4 1 5 6 7 8]) = velmat;
test_series(size(velmat,1)+1:end,1)             = loc_props(:,1);
test_series(size(velmat,1)+1:end,9:end)         = loc_props(:,2:end);

% Create matrix with dataset of properties to be interpolated (northing,
% easting, heading, velocity)
% ts_interp = test_series(isnan(test_series(:,2)),:);
ts_interp = unique(test_series(isnan(test_series(:,2)),[1 9 10 11 12]),'rows');

% For interplation, the function unique is not sufficient to delete enough
% entries. The time data will be checked for differences of 0, which
% depends on the machine precision.
delta = diff(ts_interp(:,1));
nodiff = find(delta == 0);

for k=1:length(nodiff)
    ts_interp(nodiff(k),1) = NaN;
end

ts_interp = ts_interp(~isnan(ts_interp(:,1)),:);

% Compute linear interpolation
% Criteria: If the entry of range in column 2
% is empty, it is a location point. Otherwise it is an interpolation point.
north_interp = interp1(ts_interp(:,1),ts_interp(:,2),...
    test_series(~isnan(test_series(:,2)),1));

east_interp = interp1(ts_interp(:,1),ts_interp(:,3),...
    test_series(~isnan(test_series(:,2)),1));
head_interp = interp1(ts_interp(:,1),ts_interp(:,4),...
    test_series(~isnan(test_series(:,2)),1));
vel_interp = interp1(ts_interp(:,1),ts_interp(:,5),...
    test_series(~isnan(test_series(:,2)),1));

% Save interpolated values to measurements
test_series(~isnan(test_series(:,2)),9) = north_interp;
test_series(~isnan(test_series(:,2)),10) = east_interp;
test_series(~isnan(test_series(:,2)),11) = head_interp;
test_series(~isnan(test_series(:,2)),12) = vel_interp;

% Delete all entries without a range:
test_series = test_series(~isnan(test_series(:,2)),:);

end
