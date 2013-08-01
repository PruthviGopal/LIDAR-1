function [ phase_velocity ] = phase_average2D(gridprops, range, Nphase,...
    lidar_sync, ldm_sync, save_figure, limit)
%% function phase_average2D
% [ phase_velocity ] = phase_average2D(gridprops, range, Nphase, lidar_sync, ldm_sync)
% [ phase_velocity ] = phase_average2D(gridprops, range,...
%    Nphase, lidar_sync, ldm_sync, save_figure)
% [ phase_velocity ] = phase_average2D(gridprops, range,...
%    Nphase, lidar_sync, ldm_sync, save_figure, limit)
% 
% DESCRIPTION
% The function applies a time averaging to the data selected
% during the algorithm. Based on the imported data, the user has to choose
% to time intervals which will be declared as measurement one and two. In
% addition, the averaging takes place according to the discrete phases
% defined in Nphase As output it provides a struct with velocity components
% dependend on their interpolation method. In addition, plots will be saved
% to folder rootfigures (see global variables) if save_figure is TRUE.
%
% INPUT
% - gridprops: array with properties of final grid: [xmin, xmax, ymin, ymax, stepsize] 
% - range: radial distance between 2 measurements. Is used to resolve the
% range gate of the LIDAR system
% - Nphase: No. of discrete sections in phase
% - lidar_sync: sync time in s to add to the data in LIDAR and location
% - ldm_sync: sync time in s to add to the data in LDM
% - save_figure: boolean. If TRUE, figures will be saved. Default is FALSE
% - limit: absolute value of maximal wind speed to be shown. Other values
% will be cropped and set by NaN in the final plots. None by default
%
% OUTPUT
% - phase_velocity: struct with data of the velocity components
%       phase_range: range of phases in rad
%       x_meshgrid: x-coordinates in meshgrid format
%       y_meshgrid: y-coordinates in meshgrid format
%       u_lin: horizontal wind component (linear)
%       v_lin: vertical wind component (linear)
%       velocity_2D_lin: absolute value of wind component (linear)
%       u_near: horizontal wind component (nearest neighbour)
%       v_near: vertical wind component (nearest neighbour)
%       velocity_2D_near: absolute value of wind component (nearest neighbour)
%       u_nat: horizontal wind component (natural neighbour)
%       v_nat: vertical wind component (natural neighbour)
%       velocity_2D_nat: absolute value of wind component (natural neighbour)
%
% Code by: Markus Schmidt
%
% $Revision: 0.6$ $Date: 2013/05/14 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

close all

% Global variables
interp_method = 3;          % 1: Linear, 2: Nearest N.,3: Natural N.
di_limit = 1.01;            % Limit for doppler intensity
vr_limit = 40;              % Limit of radial velocity
known_times = true;     % Times of measurements are known
figure_limit = 12;          % Limit of velocity shown in figures

% Input Check
if nargin > 7 || nargin < 5
    error('Incorrect number of input arguments')
end

if ~exist('save_figure','var')
    save_figure = false;
end

if ~exist('limit','var')
    is_limit = false;
else
    is_limit = true;
end

scriptfolder = pwd;

%% Import of data
test_series = import_datav2(range,lidar_sync,ldm_sync);

%% Select measurements
% Manually defined start and endpoint for testing the algorithm. UI is
% written below.
if known_times
    start1 = [19, 00, 00];  %6.61e4, but missing GPS data
    end1 = [19, 33, 20];    %6.68e4
    start2 = [19, 38, 20];  %6.71e4
    end2 = [20, 26, 40];    %7e4
    [ measure_1, measure_2 ] = ts_selector( test_series, start1, end1,...
        start2, end2);
else
    [ measure_1, measure_2 ] = ts_selector(test_series);
end

%% Evaluate distance between the 2 measurements
% Averaging the values for latitude and longitude, ignoring entries with
% NaN
lat1_mean = mean(measure_1(~isnan(measure_1(:,9)),9));
lon1_mean = mean(measure_1(~isnan(measure_1(:,10)),10));
lat2_mean = mean(measure_2(~isnan(measure_2(:,9)),9));
lon2_mean = mean(measure_2(~isnan(measure_2(:,10)),10));

% Compute distance between the 2 positions. Since UTM is in kilometres, a
% dimension change to metres takes place:
distance = sqrt((lat1_mean-lat2_mean)^2+(lon1_mean-lon2_mean)^2);

% Ask the user if where the 2nd measurement is situated
asking = true;
while asking
    disp('The relative position of the two measurements has to be defined.')
    prompt = 'Is measurement 2 situated in positive or negative x-direction regarding measurement 1? [P/N]:';
    orientation = input(prompt,'s');
    
    if orientation == 'n' || orientation == 'N'
        distance = -distance;
        asking = false;
    elseif orientation == 'p' || orientation == 'P'
        asking = false;
    else
        disp('Wrong input, please try again!')
    end
end

clear lat1_mean lon1_mean lat2_mean lon2_mean orientation  % clear unused variables

%% Compute phase information of measurements
[ measure_phase_1 ] = phase_finder( measure_1 );
[ measure_phase_2 ] = phase_finder( measure_2 );

%% Filter data
% Filter all data according to limit of doppler intensity
measure_phase_1 = measure_phase_1(measure_phase_1(:,4) >= di_limit,:);
measure_phase_2 = measure_phase_2(measure_phase_2(:,4) >= di_limit,:);
% Filter all ranges below 50 metres
measure_phase_1 = measure_phase_1(measure_phase_1(:,2) > 50,:);
measure_phase_2 = measure_phase_2(measure_phase_2(:,2) > 50,:);
% Apply filter for all radial wind speeds above 30 m/s
measure_phase_1 = measure_phase_1(measure_phase_1(:,3) < vr_limit & ...
    measure_phase_1(:,3) > -vr_limit,:);
measure_phase_2 = measure_phase_2(measure_phase_2(:,3) < vr_limit &...
    measure_phase_2(:,3) > -vr_limit,:);

%% Interpolation
for k=0:Nphase-1
    phase_range = [k, k+1] *2/3*pi/Nphase;
    if is_limit
    phase_velocity(k+1) = phase_interp( measure_phase_1, measure_phase_2,...
        gridprops, phase_range, distance, limit);
    else
    phase_velocity(k+1) = phase_interp( measure_phase_1, measure_phase_2,...
        gridprops, phase_range, distance);
    end 
end

%% Visualisation of results
if save_figure
    % Define location to save figures
    if exist('rootfigures.txt','file')
        import = importdata('rootfigures.txt');
        rootfigures = cell2mat(import(1,1));
    else
        % UI for selecting the data folders
        disp('Please select the root folder for figures.')
        rootfigures = uigetdir('','Please select the root folder for figures.');
        
        % Save path to txt-file
        fileID = fopen('rootfigures.txt','wt');
        fprintf(fileID, '%s\n', rootfigures);
        fclose(fileID);
    end
    
    % Check existence of folder figures. If not present, create one
    latest = 1;
    if exist(rootfigures,'dir')~=7
        mkdir(rootfigures);
    end  
    cd(rootfigures);
    list = dir('phase_wind_speed-*');
    list = struct2cell(list);
    list = list(1,:);           % Crop list to filenames
    list = char(list);
    A = NaN(size(list,1),1);
    for k = 1:size(list,1)
        temp = strsplit(list(k,:), {'-', '.'},'CollapseDelimiters',true);
        A(k) = str2num(cell2mat(temp(2)));
    end
    latest = max(A) + 1;
    if isempty(latest)
        warning('Existing plots-folder in figures could not be found.')
        latest = 1;
    end
    currentfolder = sprintf('%s-%02d','phase_wind_speed',latest);
    mkdir(currentfolder);
    disp(['Figures will be saved in folder number ', num2str(latest), '.']);
    
    % Function handels
    printeps = @(fname) print('-depsc2',... % print eps file
        fullfile(rootfigures,currentfolder,fname),'-zbuffer','-r200');
    printpng = @(fname) print('-dpng',... % print png file
        fullfile(rootfigures,currentfolder,fname),'-zbuffer','-r200');
    savefig = @(fname) hgsave(gcf,...       % Save figure
        fullfile(rootfigures,currentfolder,fname));
    
    set(0,...                               % figure settings
        'DefaultFigureColormap',gray,...
        'DefaultAxesVisible','on',...
        'DefaultAxesNextPlot','add',...
        'DefaultAxesFontSize',10);

    % 2D velocities with linear interpolation
    switch interp_method
        case 1
            titlestart = 'wind speed [m/s] in 2D plane (Linear Interpolation) for ';
            namestart = 'phase_vel_lin';
            velocity = NaN(size(phase_velocity(1).x_meshgrid,1),...
                size(phase_velocity(1).x_meshgrid,2),Nphase);
            for m=1:Nphase
                if ~isempty(phase_velocity(m).velocity_2D_lin)
                    velocity(:,:,m) = phase_velocity(m).velocity_2D_lin;
                else
                    velocity(:,:,m) = [];
                end
            end
        case 2
            titlestart = 'wind speed [m/s] in 2D plane (Nearest N.) for ';
            namestart = 'phase_vel_near';
            velocity = NaN(size(phase_velocity(1).x_meshgrid,1),...
                size(phase_velocity(1).x_meshgrid,2),Nphase);
            for m=1:Nphase
                if ~isempty(phase_velocity(m).velocity_2D_near)
                    velocity(:,:,m) = phase_velocity(m).velocity_2D_near;
                else
                    velocity(:,:,m) = [];
                end
            end
        case 3
            titlestart = 'wind speed [m/s] in 2D plane (Natural N.) for ';
            namestart = 'phase_vel_nat';
            velocity = NaN(size(phase_velocity(1).x_meshgrid,1),...
                size(phase_velocity(1).x_meshgrid,2),Nphase);
            for m=1:Nphase
                if ~isempty(phase_velocity(m).velocity_2D_nat)
                    velocity(:,:,m) = phase_velocity(m).velocity_2D_nat;
                else
                    velocity(:,:,m) = [];
                end
            end
    end
    for l=1:Nphase
        if ~isnan(velocity(:,:,l))
        hfigure(gcf) = figure(gcf+1);
        contourf(phase_velocity(l).x_meshgrid,...
            phase_velocity(l).y_meshgrid,velocity(:,:,l))
        title([titlestart,...
            num2str(phase_velocity(l).phase_range(1),3), ' to ',...
            num2str(phase_velocity(l).phase_range(2),3), 'radians.']);
        xlabel('horizontal distance to LIDAR in m')
        ylabel('vertical distance to LIDAR in m')
        if is_limit
            caxis([0 limit])
        else
            caxis([0 figure_limit])
        end
        colorbar
        colormap('jet')
        name = sprintf('%s%02d',namestart,l);
        printeps(name);
        printpng(name);
        savefig(name);
        else
            % If no data is present, create .txt file dummy to indicate it.
            name = sprintf('%s%02d%s',namestart,l,'.txt');
            fid = fopen(fullfile(rootfigures,currentfolder,name), 'w');
            fprintf(fid, '%s', 'No interpolation available');
            fclose(fid); 
        end
    end
end
close all
cd(scriptfolder)
end
