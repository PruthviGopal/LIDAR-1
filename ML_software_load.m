function [ loc_props ] = ML_software_load(start_time,end_time, log_home )
%% function location_load
% function [ loc_props ] = ML_software_load(start_time,end_time, log_home )
% 
% DESCRIPTION
% This function is used to import the GPS data aquired by the WindRover
% during measurements. In addition, the data is converted from GPS
% coordinates to UTM(Universal Transverse Mercator coordinate system).
% 
% INPUT
% - start_time: numerical array of pattern [yyyy, m, d, h, min, sec]
% - end_time:   numerical array of pattern [yyyy, m, d, h, min, sec]
% - log_home:   string containing the location folder of log files.
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
% $Revision: 0.4$ $Date: 2013/04/22 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Input check
if nargin~=3
    error('Wrong number of input arguments.')
end

% Check with operating system is used. Important to create full path names.
if isunix || ismac
    slash = '/';
else
    slash = '\';
end

% Save current folder
script_folder = pwd;

% Import of file:
cd(log_home)
filename = sprintf('%4d%s%02d%s%02d%s',start_time(1),' ',start_time(2),' ',...
    start_time(3),'*.txt');
log_listing = dir(filename);
if isempty(log_listing)
    error('Could not find .txt files in assigned folder.')
end
temp = pwd;
file_list = [];
for p=1:size(log_listing,1)
    file_list = [file_list; sprintf('%s%s%s',temp, slash,log_listing(p).name)];
end

% Compute the .txt files
cd(script_folder)
gps_load = [];
for k=1:size(file_list,1)
gps_load = [gps_load;load(file_list(k,:))];
end

% Filter data and use only entries with fix = true
gps_load = gps_load(gps_load(:,17) == true,:);

% Calculate the daytime in s
time = 3600*gps_load(:,4)+60*gps_load(:,5) + gps_load(:,6) + (10e-3)*gps_load(:,7);

% Limit the data to start and end time
time_min = 3600*start_time(:,4)+60*start_time(:,5) + start_time(:,6);
time_max = 3600*end_time(:,4)+60*end_time(:,5) + end_time(:,6);
time_limit = boolean((time>=time_min).*(time<=time_max));
gps_load = gps_load(time_limit,:);
time = time(time_limit,:);

% Column properties of gps_load and its allocation to the loc_props:
% 01 systime yyyy       -> none
% 02 systime mm         -> none
% 03 systime dd         -> none
% 04 systime hh         ->1
% 05 systime minmin     ->1
% 06 systime ss         ->1
% 07 systime msmsms     ->1
% 08 gpstime dd         -> none
% 09 gpstime hh         -> none
% 10 gpstime minmin     -> none
% 11 gpstime ss         -> none
% 12 gpstime msmsms     -> none
% 13 longitude          -> 2 (GPS2UTM)
% 14 latitude           -> 3 (GPS2UTM)
% 15 altitude           -> none
% 16 sat                -> none
% 17 fix                -> none
% 18 gpsCog             -> 4
% 19 gpssog             -> 5
% 20 avyaw              -> none
% 21 avroll             -> none
% 22 avpitch            -> none

% Convert the GPS data to UTM
northing = NaN(size(gps_load,1),1);
easting = NaN(size(gps_load,1),1);
for i=1:size(gps_load,1)
    [northing(i), easting(i)] = GPS2UTM(gps_load(i,14)*10e-8, gps_load(i,13)*10e-8);
end

% Find fragments in data. If there is a time gap, the missing of data is
% assumed. Filtering will only take place within these fragments
frag_diff = diff(time);
frag_true = frag_diff >=10*median(frag_diff);
frag_pos = find(frag_true);
frag_pos = [0; frag_pos; length(time)];

% Check if WindRover is moving. If yes, the heading can be applied. If no,
% a heading has to be found according to the mean of the fragment
ismoving = gps_load(:,19) >0;

%The heading is defined as 0 if heading north.
heading = NaN(size(northing));
heading(ismoving,1) = gps_load(ismoving,18);

inrow = [];
for i=1:length(heading)
    if ~ismoving
        inrow = [inrow; i];
    elseif length(inrow) < 10
        heading(inrow,1) = NaN;
    else
        heading(inrow,1) = median(gps_load(inrow,18));
        inrow = [];
    end
end

% % Compute running average of GPS data fragments
% windowSize = 50;
% n_ra = [];
% e_ra = [];
% valid_time = [];
% for i=1:length(frag_pos)-1
%     % For a proper fragment size, the measurement row has to be at least 10
%     % entries long
%     if (frag_pos(i+1)-frag_pos(i)+1) >= windowSize
%         n_ra_temp = filter(ones(1,windowSize)/windowSize,1,...
%             northing(frag_pos(i)+1:frag_pos(i+1),1));
%         e_ra_temp = filter(ones(1,windowSize)/windowSize,1,...
%             easting(frag_pos(i)+1:frag_pos(i+1),1));
%         n_ra = [n_ra; n_ra_temp];
%         e_ra = [e_ra; e_ra_temp];
%         valid_time(frag_pos(i)+1:frag_pos(i+1),1) = true;
%     else
%         valid_time(frag_pos(i)+1:frag_pos(i+1),1) = false;
%     end
% end
% 
% % Reduce time entries
% time_ra = time(boolean(valid_time),1);
% 
% % Filter results of running average for outyers
% n_valid = n_ra > 0.999*mean(n_ra);
% e_valid = e_ra > 0.999*mean(e_ra);
% mean_filter = boolean(n_valid.*e_valid);
% time_ra = time(mean_filter);
% n_ra = n_ra(mean_filter);
% e_ra = e_ra(mean_filter);
% 
% northing_fit = fit(time_ra,n_ra,'smoothingspline','SmoothingParam',1);
% easting_fit = fit(time_ra,e_ra,'smoothingspline','SmoothingParam',1);
% % Define times for the differentiation between the measurement points:
% delta_t = time_ra(1:end-1,1) + 0.5*diff(time_ra);
% 
% delta_n = differentiate(northing_fit,delta_t);
% delta_e = differentiate(easting_fit,delta_t);
% 
% hold on
% scatter(time_ra,n_ra,'x')
% plot(northing_fit)

% % Apply a linear regression to the location data
% G = [ones(size(time)), time];
% popt_n = G \ northing;
% popt_e = G \ easting;
% northing = G * popt_n;
% easting = G * popt_e;

% Compute the current velocity. 
velocity = NaN(size(time,1),1);
delta_t = diff(time);
delta_n = diff(northing);
delta_e = diff(easting);
velocity(1:end-1) = sqrt(delta_n.^2+ delta_e.^2)./delta_t;
velocity(end) = velocity(end-1);

% Filter the velocity for outlyers
velocity(abs(5*mean(velocity))< abs(velocity))= NaN;

% % Apply a runnning average to heading and velocity
% windowSize = 50;
% vel_ra(vel_in) = filter(ones(1,windowSize)/windowSize,1,velocity(vel_in));
% vel_ra(~vel_in) = NaN;

% Set up the final matrix loc_props
loc_props = [time, northing, easting, heading, velocity];

% Filter results to get meaningful data
loc_props = loc_props(~isnan(loc_props(:,4)),:);

end

