function [ v_r_data ] = interp_radial2D( test_series, save_figure, doppler_i )
%% function interp_radial2D
% [ v_r_data ] = interp_radial2D( test_series, save_figure, doppler_i )
% 
% DESCRIPTION
% The function compute a time average of all provides data without taking
% into acount phase dependency of values. It takes a measurement on a 
% constant position as input. Out of this measurement, it creates an
% interpolated function of class 'scatteredInterpolant' as output.
% If save_figures is provided with value TRUE, figures of several
% properties will be saved to subfolder "figures". No input of save_figures
% will result in no export of figures.
%
% The following filters are applied to the data:
% - Filtering by statistical meanings
%   Assuming a alpha= 0.05, the confidence interval is 95%. If this interval
%   yields a broader deviation of the mean than 10 %, the averaged value is
%   replaced by NaN and will be interpolated.
%   transform((0.1.*abs(transform(:,3)) >= 2.131*sqrt(transform(:,4))./sqrt(transform(:,5))),3) = NaN;
%
% INPUT
% - test_series: 
%    Column 1: Time in seconds, starting at 12AM of the measurement day
%    column 2: range in m. Radial distance from LIDAR
%    column 3: Radial velocity in m/s
%    column 4: Doppler intensity
%    column 5: azimuthal angle in degrees
%    column 6: elevation angle in degrees
%    column 7: pitch in ??
%    column 8: roll in ??
%    Column 9: Norhting in m of LIDAR
%    Column 10: Easting in m of LIDAR
%    Column 11: The computed number of revolutions (nor) in seconds.
% - save_figure: boolean. If TRUE, figures will be saved. Default is FALSE
% - doppler_i: integer doppler intensity, is used as lower threshold of
% data, default value: 1.01
% 
% OUTPUT
% - v_r_data: struct with the following entries
%       x_meshgrid: x-coordinates created by function 'meshgrid'
%       y_meshgrid: y-coordinates created by function 'meshgrid'
%       linear_interp: object of class 'scatteredInterpolant'
%       nearest_nb: object of class 'scatteredInterpolant'
%       natural_nb:  object of class 'scatteredInterpolant'
%       elevation: elevations according to grid coordinates
%
% Code by: Markus Schmidt
%
% $Revision: 1.1$ $Date: 2013/05/13 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

close all
% Check input variable
if size(test_series,2)< 11 || size(test_series,2)> 12
    error('Input does not have the right dimensions. Please try again!')
end

% Check number of input arguments
if nargin<2
    save_figure = false;
end

if nargin<3
    doppler_i = 1.01;
end

%% Compute ranges and elevations
% Find all ranges and elevation angles (in radians) present in the measurement. NaN
% number are negleted. Finally, the 'unique' function combines both
% outputs.
ranges = test_series(~isnan(test_series(:,2)),2);
el = test_series(~isnan(test_series(:,6)),6);

% Delete double entries, change to radians
ranges = unique(ranges);
el = unique(deg2rad(el));

%% Check existence of folder figures. If not present, create one
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
    
    scriptfolder = pwd;
    latest = 1;
    if exist(rootfigures,'dir')~=7
        mkdir(rootfigures);
    elseif size(dir(rootfigures),1) > 2   % Check wether folder is empty
        % Check numbering of figures. If figures are already present, the number
        % should increase by one number
        cd(rootfigures);
        list = dir('range_elevation*');
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
end

%% Create plots of range,elevation vs. doppler intensity for quality analysing
if save_figure
    hfigure(1) = figure(1);
    set(gcf,'PaperUnits', 'centimeters', 'PaperType', 'A4','PaperOrientation','portrait','PaperPositionMode', 'manual','PaperPosition', [2 1 27.7 21])
    
    % Plot range vs. doppler intensity
    subplot(1,2,1)
    scatter(test_series(:,2),test_series(:,4))
    axis([0,max(ranges),0.6,1.3])
    title(sprintf('%s%02d','Range vs. doppler intensity for measurement ',latest))
    xlabel('Range in m')
    ylabel('Doppler intensity')
    
    % Plot elevation vs. doppler intensity
    subplot(1,2,2)
    scatter(test_series(:,6),test_series(:,4))
    axis([0,rad2deg(max(el)),0.6,1.3])
    title(sprintf('%s%02d','Elevation vs. doppler intensity in ',latest))
    xlabel('Elevation in degrees')
    ylabel('Doppler intensity')
    
    name = 'range_elevation_doppler';
    printeps(name);
    savefig(name);
end

%% Filter and average data

% Filter data with 3 conditions
% 1) range should not be NaN
% 2) elevation should not be NaN
% 3) Measurement should reach a Doppler Intensity of at doppler_i
% 4) All ranges below 50 metres are neglected

filter = test_series((~isnan(test_series(:,2))),:);     % Filter range NaNs
filter = filter((~isnan(filter(:,3))),:);               % Filter elev. NaNs
filter = filter((filter(:,4)>=doppler_i),:);            % Filter bad data
filter = filter(filter(:,2) > 50,:);     % Filter range below 50 metres

% Crop result to save memory. Only range, radial wind speed and elevation
% is maintained. The elevation is changed from degree to radian
filter = filter(:,[2, 3, 6]);
filter(:,3) = deg2rad(filter(:,3));

% Calculate average values for each measurement and each coordinate
% Structure of matrix 'averg':
% Column 1: range in m
% Column 2: elevation in rad
% Column 3: mean radial wind speed
% Column 4: Standard deviation of mean wind speed
% Column 5: Number of measurements

averg = NaN(length(ranges)*length(el),5);

% Collection of all measurements of one location (same range and
% elevation).
n=1;
for k=1:length(ranges)
    for l=1:length(el)
        temp = filter((filter(:,1) == ranges(k)),:);        % Filter one range
        temp = temp((temp(:,3) == el(l)),:);            % Filter one elevation
        if isempty(temp)
            averg(n,:) = NaN;
        else
            averg(n,:) = [temp(1,1), temp(1,3), mean(temp(:,2)),...
                std(temp(:,2)), size(temp,1)]; % Save mean value of result
        end
        n = n+1;
    end
end

% Delete NaN entries in mean
averg = averg(~isnan(averg(:,1)),:);

clear filter temp   % clear unsused variables

% Apply coordinate system transformation from cylindrical to carthesian,
% store the x and y coordinate in new columns 6 and 7.
% Structure of matrix 'transform':
% Column 1: range in m
% Column 2: elevation in rad
% Column 3: mean radial wind speed
% Column 4: Standard deviation of mean wind speed
% Column 5: Number of measurements
% Column 6: x-coordinate in grid
% Column 7: y-coordinate in grid

transform = [averg, averg(:,1).*cos(averg(:,2)), averg(:,1).*sin(averg(:,2))];


%% Create Plots of the point averaged measurements as a velocity field
% Print the radial velocities of the 2 measurements
% Plot elevation vs. doppler intensity
if save_figure
    hfigure(gcf) = figure(gcf+1);
    set(gcf,'PaperUnits', 'centimeters', 'PaperType', 'A4','PaperOrientation','portrait','PaperPositionMode', 'manual','PaperPosition', [2 1 27.7 21])
    
    quiver(transform(:,6),transform(:,7),transform(:,3).*cos(transform(:,2)),transform(:,3).*sin(transform(:,2)));
    title(sprintf('%s%02d','radial velocity measurement ',latest))
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')
    
    name = 'quiver_unfiltered';
    printeps(name);
    savefig(name);
    
    % Statistical values of measurements
    hfigure(gcf) = figure(gcf+1);
    set(gcf,'PaperUnits', 'centimeters', 'PaperType', 'A4','PaperOrientation','portrait','PaperPositionMode', 'manual','PaperPosition', [2 1 27.7 21])
    
    subplot(1,2,1);
    scatter(transform(:,1),transform(:,4));
    title(sprintf('%s%02d','Standard deviation measurement ',latest))
    % axis([0,max(ranges),0, max(transform(:,4))])
    xlabel('radial distance to LIDAR in m')
    ylabel('standard deviation')
    
    subplot(1,2,2);
    scatter(transform(:,1),transform(:,5));
    title(sprintf('%s%02d','Number of valid measurements in ',latest))
    % axis([0,max(ranges),0, max(transform(:,5))])
    xlabel('horizontal distance to LIDAR in m')
    ylabel('number of measurements')
    
    name = 'statistics';
    printeps(name);
    savefig(name);
end

% Create plots with interpolated values of standard deviation
if save_figure
    std_grid = griddata(transform(:,6), transform(:,7), transform(:,4),...
        min(transform(:,6)):max(transform(:,6)),...
        (min(transform(:,7)):max(transform(:,7)))');
    
    hfigure(gcf) = figure(gcf+1);
    set(gcf,'PaperUnits', 'centimeters', 'PaperType', 'A4','PaperOrientation','portrait','PaperPositionMode', 'manual','PaperPosition', [2 1 27.7 21])
    
    contourf(min(transform(:,6)):max(transform(:,6)),(min(transform(:,7)):max(transform(:,7)))',std_grid);
    title(sprintf('%s%02d','Standard deviation (Linear Interpolation) in measurement ',latest))
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')
    colorbar('location','southoutside')
    colormap jet
    
    name = 'standard deviation';
    printeps(name);
    savefig(name);
end

%% Filtering by statistical meanings
% % Assuming a alpha= 0.05, the confidence interval is 95%. If this interval
% % yields a broader deviation of the mean than 10 %, the averaged value is
% % replaced by NaN and will be interpolated.
% transform((0.1.*abs(transform(:,3)) >= 2.131*sqrt(transform(:,4))./sqrt(transform(:,5))),3) = NaN;
% 
% %% Plot the results of statistical filtering
% if save_figure
%     hfigure(gcf) = figure(gcf+1);
%     set(gcf,'PaperUnits', 'centimeters', 'PaperType', 'A4','PaperOrientation','portrait','PaperPositionMode', 'manual','PaperPosition', [2 1 27.7 21])
%     quiver(transform(:,6),transform(:,7),transform(:,3).*cos(transform(:,2)),transform(:,3).*sin(transform(:,2)));
%     hold on
%     plot((transform(isnan(transform(:,3)),6)),(transform(isnan(transform(:,3)),7)),'.','Color','r','MarkerSize',6);
%     hold off
%     
%     title(sprintf('%s%02d%s','radial v measurement ',latest,'(\alpha = 0.05) '))
%     xlabel('horizontal distance to LIDAR in m')
%     ylabel('vertical distance to LIDAR in m')
%     
%     name = 'quiver_filtered';
%     printeps(name);
%     savefig(name);
% end

%% Interpolation
% Use the scatteredInterpolant class to set up interpolated data. The field
% is limited to the defined x and y values to avoid extrapolation of
% unnecessary areas
x_min = min(transform(:,6));
x_max = max(transform(:,6));
y_min = min(transform(:,7));
y_max = max(transform(:,7));

% Filter NaN created by statistical filter
transform = transform(~isnan(transform(:,3)),:);

% Linear Interpolation
lin_interp = scatteredInterpolant(transform(:,6),transform(:,7),transform(:,3));

% Nearest Neighbour interpolation
near_interp = scatteredInterpolant(transform(:,6),transform(:,7),transform(:,3),'nearest');

% Natural Neighbour interpolation
nat_interp = scatteredInterpolant(transform(:,6),transform(:,7),transform(:,3),'natural');

% Create new data matrix and save entries according to grid
% Structure of matrix 'griddata_*': The matrix contains an interpolated
% value for each 
[xfinal_meshgrid, yfinal_meshgrid] = meshgrid(x_min:x_max, y_min:y_max);

griddata_lin = lin_interp(xfinal_meshgrid, yfinal_meshgrid);
griddata_near = near_interp(xfinal_meshgrid, yfinal_meshgrid);
griddata_nat = nat_interp(xfinal_meshgrid, yfinal_meshgrid);

%% Calculation of elevation grid
el_grid = atan2(yfinal_meshgrid,xfinal_meshgrid);

%% Create output struct
v_r_data = struct('x_meshgrid',xfinal_meshgrid,...
    'y_meshgrid',yfinal_meshgrid,'linear_interp',lin_interp,...
    'nearest_nb', near_interp,'natural_nb', nat_interp,...
    'elevation', el_grid);

%% Plot linear interpolated radial wind speeds
if save_figure
    hfigure(gcf) = figure(gcf+1);
    set(gcf,'PaperUnits', 'centimeters', 'PaperType', 'A4','PaperOrientation','portrait','PaperPositionMode', 'manual','PaperPosition', [2 1 27.7 21])
    
    contourf(xfinal_meshgrid, yfinal_meshgrid, griddata_lin);
    title('Radial wind speed (Linear Interpolation) measurement 1 ')
    title(sprintf('%s%02d','Radial wind speed (Linear Interpolation) measurement ',latest))
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')
    colorbar
    colormap('jet')
    
    name = 'v_r_lin_interp';
    printeps(name);
    savefig(name);
    
    % Plot nearest neighbor interpolated radial wind speeds
    hfigure(gcf) = figure(gcf+1);
    set(gcf,'PaperUnits', 'centimeters', 'PaperType', 'A4','PaperOrientation','portrait','PaperPositionMode', 'manual','PaperPosition', [2 1 27.7 21])
    
    contourf(xfinal_meshgrid, yfinal_meshgrid,griddata_near);
    title(sprintf('%s%02d','Radial wind speed (Nearest Neighbour) measurement ',latest))
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')
    colorbar
    colormap('jet')
    
    name = 'v_r_near_interp';
    printeps(name);
    savefig(name);
    
    % Plot natural neighbor interpolated radial wind speeds
    hfigure(gcf) = figure(gcf+1);
    set(gcf,'PaperUnits', 'centimeters', 'PaperType', 'A4','PaperOrientation','portrait','PaperPositionMode', 'manual','PaperPosition', [2 1 27.7 21])
    
    contourf(xfinal_meshgrid, yfinal_meshgrid,griddata_nat);
    title(sprintf('%s%02d','Radial wind speed (Natural Neighbour) measurement ',latest))
    xlabel('horizontal distance to LIDAR in m')
    ylabel('vertical distance to LIDAR in m')
    colorbar
    colormap('jet')
    
    name = 'v_r_nat_interp';
    printeps(name);
    savefig(name);
end

close all
end
