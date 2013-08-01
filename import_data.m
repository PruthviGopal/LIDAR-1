function [ test_series ] = import_data( range )
%% function import_data
% function [ test_series ] = import_data( range )
% 
% DESCRIPTION The function creates a n times 11 matrix out of the provides
% measurement files. For LIDAR data, only the main folder has to be
% supplied. LDM and GPS data has to be provided in .txt files.
%
% INPUT
% - range: radial distance between 2 measurements. Is used to resolve the
% range gate of the LIDAR system
%
% OUTPUT
% - test_series: matrix with the following data
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
% Code by: Markus Schmidt
%
% $Revision: 1.1$ $Date: 2013/04/25$
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Error messages
if nargin ~= 1
    error('Incorrect number of input arguments')
end

% Save current folder
script_folder = pwd;

%% Selecting the data folders
if exist('lastpath.txt','file')
    import = importdata('lastpath.txt');
    lidar_folder = cell2mat(import(1,1));
    ldm_file = cell2mat(import(2,1));
    loc_file = cell2mat(import(3,1));
else
    % UI for selecting the data folders
    disp('Please select the main folder with LIDAR data.')
    lidar_folder = uigetdir('','Please select the day folder with LIDAR data:');
    
    disp('Please select the file with LDM data.')
    [ldmName,ldmPath,~] = uigetfile('*.txt','Please select the file with LDM data:');
    ldm_file = fullfile(ldmPath, ldmName);
    
    disp('Please select the file with GPS data.')
    [locName,locPath,~] = uigetfile('*.txt','Please select the file with GPS data:');
    loc_file = fullfile(locPath, locName);
    
    % Save path to txt-file
    fileID = fopen('lastpath.txt','wt');
    fprintf(fileID, '%s\n', lidar_folder);
    fprintf(fileID, '%s\n', ldm_file);
    fprintf(fileID, '%s\n', loc_file);
    fclose(fileID);
end

%% Compute the available data

% With the given main folder of LIDAR data, the algorithm goes into each
% subfolder and extracts the given scn files
cd(lidar_folder)
velmat = [];
folder_listing = dir;

% Create a list with all file locations
file_list = [];
for k=3:size(folder_listing,1)
    if folder_listing(k).isdir
        cd(folder_listing(k).name)
        lidar_listing = dir('*.scn');
        temp = pwd;
        for m=1:size(lidar_listing,1)
            file_list = [file_list; sprintf('%s%s%s',temp, slash,lidar_listing(m).name)];
        end
    end
    cd('../')
end

% Compute the listed files with vel_load
cd(script_folder)
for m=1:size(file_list,1)
    velmat = [velmat; vel_load(file_list(m,:), range)];
end

% Load the other data
location = location_load(loc_file);
nor = blade_nor(ldm_file);

test_series = fuse_data(velmat, location, nor);
end

