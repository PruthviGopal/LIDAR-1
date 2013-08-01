function [ measure_phase ] = phase_finder( measure )
%% function phase_finder
% function [ measure_phase ] = phase_finder( measure )
% 
% DESCRIPTION
% The function computes the phase of a the LIDAR measurements present in
% the input measure. Assuming a discretized number of revolutions, by the
% time shift between two data point the phase is computed. The computed
% phase in radians is added as a 12th column to the dataset.
%
% INPUT
% - measure: 
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
%    Column 11: The computed number of revolutions (nor) in seconds.
%
% OUTPUT
% - measure_phase: 
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
%    Column 11: The computed number of revolutions (nor) in seconds.
%    Column 12: Phase of blade in radians
%
% Code by: Markus Schmidt
%
% $Revision: 1.1$ $Date: 2013/05/13$
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Check input variable
if size(measure,2)< 11 || size(measure,2) > 12 
    error('Input of measurements do not have the right dimensions. Please try again!')
end

if nargin ~= 1
    error('You must provide 1 input argument.')
end

% Divide the input into LIDAR and nor information
nor = measure(~isnan(measure(:,11)),:);
speed = measure(~isnan(measure(:,2)),:);

if isempty(nor) || isempty(speed)
    error(' Input measurement does not contain any LIDAR data.')
end

% Compute the phase information for each LIDAR measurement
phase = NaN(length(speed),1);
for k=1:length(speed)
    rrow = find(nor(:,1) < speed(k,1));
    if isempty(rrow)
        phase(k,1) = NaN;
    else
        phase(k,1) = (speed(k,1)-nor(rrow(end,1),1))*nor(rrow(end,1),11);
    end
end

% Filter data
% Values have to be below 2/3*pi, since the wind turbine has 3 blades which
% result in this maximal angle
speed_c = [measure(~isnan(measure(:,2)),:), phase];
speed_c = speed_c(phase<=2/3*pi,:);

% Create output variable
measure_phase = [speed_c; nor, NaN(length(nor),1)];
measure_phase = unique(measure_phase,'rows');

end

