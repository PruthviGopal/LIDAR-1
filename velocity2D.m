function [ vel_comps ] = velocity2D( measure_1, measure_2, gridprops, save_figure, limit )
%% function velocity2D
% function [ vel_comps ] = velocity2D( measure_1, measure_2, gridprops, save_figure, limit )
% 
% DESCRIPTION
% The function takes to measurements (measure_+ and measure_2) and provides
% velocities in a 2D plane. First, the distance between the measurements is
% calculated. The assumption is, that the wake is situated between the two
% measurements, while azimuth=0° of measurement 1 points to location of
% measurement 2. Figures can be saved in subfolder figure if desired. An
% additional filter for the final velocity can be applied via input
% 'limit'.
%
% INPUT
% - measure_1, measure_2: 
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
% - gridprops: array with properties of final grid: [xmin, xmax, ymin, ymax, stepsize] 
% - save_figure: boolean. If TRUE, figures will be saved. Default is FALSE
% - limit: absolute value of maximal wind speed to be shown. Other values
% will be cropped and set by NaN in the final plots. Zero by default
%
% OUTPUT
% - vel_comps: struct with data of the velocity components
%   x_meshgrid: x-coordinates in meshgrid format
%   y_meshgrid: y-coordinates in meshgrid format
%   u_lin: horizontal wind component (linear)
%   v_lin: vertical wind component (linear)
%   u_near: horizontal wind component (nearest neighbour)
%   v_near: vertical wind component (nearest neighbour)
%   u_nat: horizontal wind component (natural neighbour)
%   v_nat: vertical wind component (natural neighbour)
%
% Code by: Markus Schmidt
%
% $Revision: 1.0$ $Date: 2013/04/05$
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Global settings
close all                               % Close all open figures

% Check input variable
if size(measure_1,2)~= 11 || size(measure_2,2)~= 11
    error('Input of measurements do not have the right dimensions. Please try again!')
end

if nargin <3
    error('You must at least provide 3 input arguments.')
end

if numel(gridprops) ~= 5
    error('Dimension of gridprops is wrong.')
end

if ~exist('save_figure','var')
    save_figure = false;
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
    prompt = 'Is measurement 2 situated in positive or negative x-direction regarding measurement 1? [p/n]:';
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

%% Interpolation
% Use function interp_radial2D for interpolation
plane1 = interp_radial2D(measure_1,save_figure);
plane2 = interp_radial2D(measure_2,save_figure);

% Create desired grid for output
xmin = gridprops(1);
xmax = gridprops(2);
ymin = gridprops(3);
ymax = gridprops(4);
stepsize = gridprops(5);

[x_mesh, y_mesh] = meshgrid(xmin:stepsize:xmax, ymin:stepsize:ymax);

% Create values from scatteredInterpolant class. Distance is applied to
% measurement to to take the shift into account
vr_lin1 = plane1.linear_interp(x_mesh, y_mesh);
vr_lin2 = plane2.linear_interp(x_mesh-distance, y_mesh);

vr_near1 = plane1.nearest_nb(x_mesh, y_mesh);
vr_near2 = plane2.nearest_nb(x_mesh-distance, y_mesh);

vr_nat1 = plane1.natural_nb(x_mesh, y_mesh);
vr_nat2 = plane2.natural_nb(x_mesh-distance, y_mesh);

% Calculate the horizontal and vertical wind speeds
u_lin = vr_lin1.*cos(atan2(y_mesh,x_mesh)) + vr_lin2...
    .*cos(atan2(y_mesh,(x_mesh-distance)));
v_lin = vr_lin1.*sin(atan2(y_mesh,x_mesh)) + vr_lin2...
    .*sin(atan2(y_mesh,(x_mesh-distance)));

u_nat = vr_nat1.*cos(atan2(y_mesh,x_mesh)) + vr_nat2...
    .*cos(atan2(y_mesh,(x_mesh-distance)));
v_nat = vr_nat1.*sin(atan2(y_mesh,x_mesh)) + vr_nat2...
    .*sin(atan2(y_mesh,(x_mesh-distance)));

u_near = vr_near1.*cos(atan2(y_mesh,x_mesh)) + vr_near2...
    .*cos(atan2(y_mesh,(x_mesh-distance)));
v_near = vr_near1.*sin(atan2(y_mesh,x_mesh)) + vr_near2...
    .*sin(atan2(y_mesh,(x_mesh-distance)));

% Compute the magnitude of velocity
velocity_2D_lin = sqrt(u_lin.^2 + v_lin.^2);
velocity_2D_near = sqrt(u_near.^2 + v_near.^2);
velocity_2D_nat = sqrt(u_nat.^2 + v_nat.^2);

% Limit the maximum and minimum of the output
if exist('limit','var')
    v_max = abs(limit);
    v_min = -abs(limit);
    
    u_lin(velocity_2D_lin > v_max) = NaN;
    v_lin(velocity_2D_lin > v_max) = NaN;
    velocity_2D_lin(velocity_2D_lin > v_max) = NaN;
    u_lin(velocity_2D_lin < v_min) = NaN;
    v_lin(velocity_2D_lin < v_min) = NaN;
    velocity_2D_lin(velocity_2D_lin < v_min) = NaN;
    
    u_near(velocity_2D_near > v_max) = NaN;
    v_near(velocity_2D_near > v_max) = NaN;
    velocity_2D_near(velocity_2D_near > v_max) = NaN;
    u_near(velocity_2D_near < v_min) = NaN;
    v_near(velocity_2D_near < v_min) = NaN;
    velocity_2D_near(velocity_2D_near < v_min) = NaN;
    
    u_nat(velocity_2D_nat > v_max) = NaN;
    v_nat(velocity_2D_nat > v_max) = NaN;
    velocity_2D_nat(velocity_2D_nat > v_max) = NaN;
    u_nat(velocity_2D_nat < v_min) = NaN;
    v_nat(velocity_2D_nat < v_min) = NaN;
    velocity_2D_nat(velocity_2D_nat < v_min) = NaN;
end

%% Visualisation of results
% 2D velocities with linear interpolation
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
    scriptfolder = pwd;
    latest = 1;
    if exist(rootfigures,'dir')~=7
        mkdir(rootfigures);
    elseif size(dir(rootfigures),1) > 2   % Check wether folder is empty
        % Check numbering of figures. If figures are already present, the number
        % should increase by one number
        cd(rootfigures);
        list = dir('2D_wind_*');
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
        cd(scriptfolder)
    end
    disp(['Figures will be saved with number ', num2str(latest), '.']);
    % Function handels
    printeps = @(fname) print('-depsc2',... % print eps file
        fullfile(rootfigures,sprintf('%s-%02d',fname,latest)), '-r300');
    
    savefig = @(fname) hgsave(gcf,...       % Save figure
        fullfile(rootfigures,sprintf('%s-%02d',fname,latest)));
    
    set(0,...                               % figure settings
        'DefaultFigureColormap',gray,...
        'DefaultAxesVisible','on',...
        'DefaultAxesNextPlot','add',...
        'DefaultAxesFontSize',10);

    hfigure(gcf) = figure(gcf+1);
    set(gcf,'PaperUnits', 'centimeters', 'PaperType', 'A4','PaperOrientation','portrait','PaperPositionMode', 'manual','PaperPosition', [2 1 27.7 21])
    
    subplot(1,2,1)
    contourf(x_mesh, y_mesh, velocity_2D_lin)
    title('wind speed in 2D plane (Linear Interpolation)')
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')
    colorbar
    colormap('jet')
    
    subplot(1,2,2)
    quiver(x_mesh,y_mesh,u_lin,v_lin)
    % hold on
    % plot(isnan(velocity_2D_lin),'.','Color','r','MarkerSize',6);
    % hold off
    title('wind speed in 2D plane (Linear Interpolation)')
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')
    
    name = '2D_wind_speed_lin';
    printeps(name);
    savefig(name);
    
    % 2D velocities with nearest neighbour
    hfigure(gcf) = figure(gcf+1);
    set(gcf,'PaperUnits', 'centimeters', 'PaperType', 'A4','PaperOrientation','portrait','PaperPositionMode', 'manual','PaperPosition', [2 1 27.7 21])
    
    subplot(1,2,1)
    contourf(x_mesh, y_mesh, velocity_2D_near)
    title('wind speed in 2D plane (nearest neighbor)')
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')
    colorbar
    colormap('jet')
    
    subplot(1,2,2)
    quiver(x_mesh,y_mesh,u_near,v_near)
    title('wind speed in 2D plane (nearest neighbor)')
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')

    name = '2D_wind_speed_near';
    printeps(name);
    savefig(name);
    
    % 2D velocities with natural neighbour
    hfigure(gcf) = figure(gcf+1);
    set(gcf,'PaperUnits', 'centimeters', 'PaperType', 'A4','PaperOrientation','portrait','PaperPositionMode', 'manual','PaperPosition', [2 1 27.7 21])
    
    subplot(1,2,1)
    contourf(x_mesh, y_mesh, velocity_2D_nat)
    title('wind speed in 2D plane (natural neighbor)')
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')
    colorbar
    colormap('jet')
    
    subplot(1,2,2)
    quiver(x_mesh,y_mesh,u_nat,v_nat)
    title('wind speed in 2D plane (natural neighbor)')
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')
    
    name = '2D_wind_speed_nat';
    printeps(name);
    savefig(name);
    
    % 2D Contour plot of horizontal velocity (Linear interpolation)
    hfigure(gcf) = figure(gcf+1);
    set(gcf,'PaperUnits', 'centimeters', 'PaperType', 'A4','PaperOrientation','portrait','PaperPositionMode', 'manual','PaperPosition', [2 1 27.7 21])
    
    subplot(1,2,1)
    contourf(x_mesh, y_mesh, u_lin)
    title('horizontal wind speed (linear interpolation)')
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')
    colorbar
    colormap('jet')
    
    subplot(1,2,2)
    contourf(x_mesh,y_mesh,v_lin)
    title('vertical wind speed (linear interpolation)')
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')
    colorbar
    colormap('jet')

    name = 'component_wind_speed_lin';
    printeps(name);
    savefig(name);
    
    % 2D Contour plot of horizontal velocity (nearest neighbour)
    hfigure(gcf) = figure(gcf+1);
    set(gcf,'PaperUnits', 'centimeters', 'PaperType', 'A4','PaperOrientation','portrait','PaperPositionMode', 'manual','PaperPosition', [2 1 27.7 21])
    
    subplot(1,2,1)
    contourf(x_mesh, y_mesh, u_near)
    title('horizontal wind speed (nearest neighbor)')
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')
    colorbar
    colormap('jet')
    
    subplot(1,2,2)
    contourf(x_mesh,y_mesh,v_near)
    title('vertical wind speed (nearest neighbor)')
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')
    colorbar
    colormap('jet')
    
    name = 'component_wind_speed_near';
    printeps(name);
    savefig(name);
    
    % 2D Contour plot of horizontal velocity (natural neighbour)
    hfigure(gcf) = figure(gcf+1);
    set(gcf,'PaperUnits', 'centimeters', 'PaperType', 'A4','PaperOrientation','portrait','PaperPositionMode', 'manual','PaperPosition', [2 1 27.7 21])
    
    subplot(1,2,1)
    contourf(x_mesh, y_mesh, u_nat)
    title('horizontal wind speed (natural neighbor)')
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')
    colorbar
    colormap('jet')
    
    subplot(1,2,2)
    contourf(x_mesh, y_mesh,v_nat)
    title('vertical wind speed (natural neighbor)')
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')
    colorbar
    colormap('jet')
    
    name = 'component_wind_speed_nat';
    printeps(name);
    savefig(name);
end

%% Create output struct
vel_comps = struct('x_meshgrid',x_mesh,'y_meshgrid',y_mesh,...
    'u_lin',u_lin, 'v_lin',v_lin,...
    'u_near',u_near, 'v_near',v_near,...
    'u_nat',u_lin, 'v_nat',v_nat);

close all
end
