function [ locationmat ] = location_load(gps_data)
%% function location_load
% function [ locationmat ] = location_load( gps_data )
% 
% DESCRIPTION
% This function is used to import the GPS data aquired by the WindRover
% during measurements. In addition, the data is converted from GPS
% coordinates to UTM(Universal Transverse Mercator coordinate system).
% 
% INPUT
% - gps_data: The File in .txt format produced by the GPS device.
% 
% OUTPUT
% - locationmat: numerical array of precision double
%   Column 1: Time in seconds, starting at 12AM of the day.
%   Column 2: Northing in m of LIDAR
%   Column 3: Easting in m of LIDAR
%
% Code by: Markus Schmidt
%
% $Revision: 1.0$ $Date: 2013/03/22 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Input check
if nargin~=1
    error('Wrong number of input arguments.')
end

% Import of file:
% Independend of number of input parametres, the .scn file is imported
gpsid=fopen(gps_data);    % Open the 'gps_data' for read access and returns integer file identifier 'gpsid' 

% Import the data available in the .txt file. Precision is double (parametre %f).
locationmat=textscan(gpsid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f');

% Column properties of cell array and its allocation to the numerical
% array:
% 01 year       ->none
% 02 month      ->none
% 03 day        ->none
% 04 hour       ->1
% 05 minute     ->1
% 06 seconds    ->1
% 07 latitude   ->2
% 08 longitude  ->3
% 09 to 14 position of truck (not important by now)

% Error Check: Data must be from same day
daycheck1 = sum(locationmat{3}(1:end/2,1));
daycheck2 = sum(locationmat{3}(end/2+1:end,1));
if daycheck1 ~= daycheck2
    disp 'The measurement is recorded on different days. Please check the input file!'
end

% Convert time into seconds and save in one cell. The range is calculated
locationmat={[locationmat{4}*3600+locationmat{5}*60+locationmat{6} locationmat{7} locationmat{8}]};
locationmat=cell2mat(locationmat);    % Convert vel_mat from cell array to numeric array

% Convert the GPS data to UTM
for i=1:length(locationmat(:,1))
    [locationmat(i,2), locationmat(i,3)] = GPS2UTM(locationmat(i,2), locationmat(i,3));
end
end

