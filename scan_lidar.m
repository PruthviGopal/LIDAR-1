function [ velmat ] = scan_lidar(start_time, end_time, lidar_home, range)
%% function scan_lidar
% [ lidar_data ] = scan_lidar(start_time, end_time, lidar_home)
% 
% DESCRIPTION
% The function imports the existing data in the LIDAR home folder according
% to the given start and end time. The return is a matrix with the relevant
% data.
%
% INPUT
% - start_time: numerical array of pattern [yyyy, m, d, h, min, sec]
% - end_time:   numerical array of pattern [yyyy, m, d, h, min, sec]
% - lidar_home: string containing the path to main folder with LIDAR data
% - range: radial distance between 2 measurements. Is used to resolve the
% range gate of the LIDAR system
%
% OUTPUT
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
%
% Code by: Markus Schmidt
%
% $Revision: 1.0$ $Date: 2013/04/24$
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Input check
if nargin ~= 4
    error('Incorrect number of input arguments')
end

if ~ischar(lidar_home)
    error('lidar_home is not a string.')
end

% Check with operating system is used. Important to create full path names.
if isunix || ismac
    slash = '/';
else
    slash = '\';
end

% Selecting the data folders
lidar_folder = lidar_home;

% Save current folder
script_folder = pwd;

%% Recursion of getting a list with the file locations
% With the given main folder of LIDAR data, the algorithm goes into each
% subfolder and extracts the given scn files
cd(lidar_folder)
velmat = [];
file_list = [];

% Find entries for year
year_list = dir;
year_start = num2str(start_time(1));
year_end = num2str(end_time(1));
[y_start, y_end] = find_folders(year_start, year_end, year_list);

for k=y_start:y_end
    cd(year_list(k).name)
    % Find entries for month
    month_list = dir;
    month_start = sprintf('%d%02d',start_time(1),start_time(2));
    month_end = sprintf('%d%02d',end_time(1),end_time(2));
    [m_start, m_end] = find_folders(month_start, month_end, month_list);
    
    for l=m_start:m_end
        cd(month_list(l).name)
        % Find entries for day
        day_list = dir;
        day_start = sprintf('%d%02d%02d',start_time(1),start_time(2),start_time(3));
        day_end = sprintf('%d%02d%02d',end_time(1),end_time(2),end_time(3));
        [d_start, d_end] = find_folders(day_start, day_end, day_list);
        for m=d_start:d_end
            cd(day_list(m).name)
            % Find entries for hour
            hh_list = dir;
            hh_start = sprintf('%02d',start_time(4));
            hh_end = sprintf('%02d',end_time(4));
            [hour_start, hour_end] = find_folders(hh_start, hh_end, hh_list);
            % Save *.scn files
            for n=hour_start:hour_end
                cd(hh_list(n).name)
                lidar_listing = dir('*.scn');
                temp = pwd;
                for p=1:size(lidar_listing,1)
                    file_list = [file_list; sprintf('%s%s%s',temp,...
                        slash,lidar_listing(p).name)];
                end
                cd ..
            end
            cd ..
        end
        cd ..
    end
    cd ..
end

%% Compute the listed files with vel_load
cd(script_folder)
for m=1:size(file_list,1)
    velmat = [velmat; vel_load(file_list(m,:), range)];
end

% Reduce the entries to the according daytime
start_sec = start_time(4)*3600 + start_time(5)*60 + start_time(6);
end_sec = end_time(4)*3600 + end_time(5)*60 + end_time(6);

velmat = velmat((velmat(:,4) >= start_sec),:);
velmat = velmat(velmat(:,4) <= end_sec,:);

end

function [start, last] = find_folders(start_name, last_name, list)
for j=1:(size(list,1))
    if strcmp(list(j).name,start_name) && list(j).isdir
        start = j;
        break
    end
end
for j=1:(size(list,1))
    if strcmp(list(j).name,last_name) && list(j).isdir
        last = j;
        break
    end
end
if ~exist('start','var') || ~exist('last','var')
    error('Could not find folders.')
end
end


