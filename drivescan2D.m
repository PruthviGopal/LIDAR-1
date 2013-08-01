function [ ] = drivescan2D( start_time, end_time, lidar_home, log_home,...
    range, map_file, long, lat,scale, color_code, usegradient, limit)
%% function drivescan2D
% [ ] = drivescan2D( start_time, end_time, lidar_home, log_home,...
%    range, map_file, long, lat,scale, color_code)
% [ ] = drivescan2D( start_time, end_time, lidar_home, log_home,...
%    range, map_file, long, lat,scale, color_code, usegradient, limit)
%
% DESCRIPTION
% Function reads existing measurement files of LIDAR and GPS. It computes a
% scan image onto the map and saves the result in the folder specified via
% GUI or in root_superposer.txt.
%
% INPUT
% - start_time: numerical array of pattern [yyyy, m, d, h, min, sec]
% - end_time:   numerical array of pattern [yyyy, m, d, h, min, sec]
% - lidar_home: string containing the path to main folder with LIDAR data
% - log_home: string containing the location folder of log files.
% - range: radial distance between 2 measurements. Is used to resolve the
% range gate of the LIDAR system
% - map_file:   location of image map
% - lat:        latitude of reference point on map
% - lon:        longitude of reference point of map
% - scale:      scale as shown on map_file
% - color_code: location of .txt-file with color code
% - usegradient: boolean. If true, computation of gradient instead of
% velocity. False by default.
% - limit: limit of colorbar to be shown on plot
%
% OUTPUT
% - none
%
% Code by: Markus Schmidt
%
% $Revision: 2.1$ $Date: 2013/05/16 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Input checks
if nargin < 10 || nargin == 11 || nargin > 12 
    error('Wrong number of input arguments.')
end

% Import lidar files
% Scan and read all files within time. Returns matrix with values
[ lidar_data ] = scan_lidar(start_time, end_time, lidar_home, range);

% Import GPS data
[ loc_props ] = ML_software_loadv2(log_home, start_time, end_time );

% Interpolate the position for the time entries of the beams
[ test_series ] = UTM_interp(lidar_data, loc_props);

% Import of map
% map_calibration
[map, UTM_map_x, UTM_map_y, ~] = map_calibration_v2(map_file,lat,long,...
    scale);

% Compute the opacity on the map and create figures
if exist('limit', 'var')
    superposer_v2(test_series, map, UTM_map_x, UTM_map_y, color_code,...
        usegradient, limit);
else
    superposer_v2(test_series, map, UTM_map_x, UTM_map_y, color_code);
end

end

