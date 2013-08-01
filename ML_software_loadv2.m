function [ loc_props ] = ML_software_loadv2(log_home, start_time, end_time)
%% function location_load
% [ loc_props ] = ML_software_loadv2(log_home)
% [ loc_props ] = ML_software_loadv2(log_home, start_time, end_time);
% DESCRIPTION This function is used to import the GPS data aquired by the
% WindRover during measurements. It requires the new data format of ML
% Software Load is introduced in May 2013. In addition, the data is
% converted from GPS coordinates to UTM (Universal Transverse Mercator
% coordinate system).
% 
% INPUT
% - log_home:   string containing the adress of folder of log files.
% - start_time: numerical array of pattern [yyyy, m, d, h, min, sec]
% - end_time:   numerical array of pattern [yyyy, m, d, h, min, sec]
% 
% OUTPUT
% - loc_props: numerical array of precision double
%   Column 1: Time in seconds, starting at 12AM of the day.
%   Column 2: Northing in m of LIDAR
%   Column 3: Easting in m of LIDAR
%   Column 4: heading in rad of LIDAR (eastwards defines 0 rad)
%   Column 5: Velocity in m/s of LIDAR
%
% Code by: Markus Schmidt
%
% $Revision: 2.1$ $Date: 2013/05/15 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Input check
if nargin < 1 || nargin > 3
    error('Wrong number of input arguments.')
end

% Save current folder
script_folder = pwd;

% Import of file:
cd(log_home)
log_listing = dir('*.txt');

if isempty(log_listing)
    error('Could not find .txt files in assigned folder.')
end

temp = pwd;
file_list = [];
for p=1:size(log_listing,1)
    file_list = [file_list; fullfile(temp, log_listing(p).name)];
end

% Compute the .txt files
cd(script_folder)
gps_load = [];
for k=1:size(file_list,1)
gps_load = [gps_load;load(file_list(k,:))];
end

% Filter data and use only entries with fix = true
gps_load = gps_load(gps_load(:,10) == true,:);

% Calculate the daytime in s
time = 3600*gps_load(:,4)+60*gps_load(:,5) + gps_load(:,6);

% Column properties of gps_load and its allocation to the loc_props:
% 01 systime yyyy       -> none
% 02 systime mm         -> none
% 03 systime dd         -> none
% 04 systime hh         -> 1
% 05 systime minmin     -> 1
% 06 systime ss         -> 1
% 07 latitude           -> 2 (GPS2UTM)
% 08 longitude          -> 3 (GPS2UTM)
% 09 altitude           -> none
% 10 fix                -> none
% 11 cog                -> 4
% 12 speed              -> 5
% 13 roll               -> none
% 14 pitch              -> none

% Convert the GPS data to UTM
northing = NaN(size(gps_load,1),1);
easting = NaN(size(gps_load,1),1);
for i=1:size(gps_load,1)
    [northing(i), easting(i)] = GPS2UTM(gps_load(i,7), gps_load(i,8));
end

% Find fragments in data. If there is a time gap, the missing of data is
% assumed. Filtering will only take place within these fragments
frag_diff = diff(time);
frag_true = frag_diff >=10*median(frag_diff);
frag_pos = find(frag_true);
frag_pos = [0; frag_pos; length(time)];

% Check if WindRover is moving. If yes, the heading can be applied. If no,
% a heading has to be found according to the mean of the fragment
ismoving = gps_load(:,12) >0;

%The heading is defined as 0 if heading north.
% heading = NaN(size(northing));
% heading(ismoving,1) = gps_load(ismoving,11);
heading = NaN(size(gps_load,1),1);
heading(1,1) = gps_load(1,11);
for k=2:size(gps_load,1)
    if gps_load(k,12) ~= 0
        heading(k,1) = gps_load(k,11);
    else
       heading(k,1) = heading(k-1,1); 
    end
end

% inrow = [];
% for i=1:length(heading)
%     if ~ismoving
%         inrow = [inrow; i];
%     elseif length(inrow) < 10
%         heading(inrow,1) = NaN;
%     else
%         heading(inrow,1) = median(gps_load(inrow,11));
%         inrow = [];
%     end
% end

% Compute the current velocity. 
velocity = gps_load(:,12);

% Filter the velocity for outlyers
velocity(abs(5*mean(velocity))< abs(velocity))= NaN;

% Set up the final matrix loc_props
loc_props = [time, northing, easting, heading, velocity];

% Crop data according to start_time and end_time
if exist('start_time','var')
    start_in_s =  3600*start_time(4) + 60*start_time(5) + start_time(6);
    end_in_s =  3600*end_time(4) + 60*end_time(5) + end_time(6);
    loc_props = loc_props(loc_props(:,1) >= start_in_s & ...
        loc_props(:,1) <= end_in_s,:);
end

end

