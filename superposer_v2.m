function [ ] = superposer_v2(test_series, map, UTM_map_x, UTM_map_y,...
    color_code, usegradient, limit)
%% function superposer_v2
% [ ] = superposer_v2(test_local, map, UTM_map_x, UTM_map_y, color_code)
% [ ] = superposer_v2(test_local, map, UTM_map_x, UTM_map_y, color_code,...
%           usegradient, limit)
%
% DESCRIPTION The function takes a given set of measurements (test_local)
% and the properties of a map (map,UTM_map_x, UTM_map_y). With this
% information, it creates a set of figures, which are a combination of the
% map and one azimuthal measurement with a defined opacity and the colours
% defined by maximal and minimal speed. The figures will be saved in the
% subfolder superposed.
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
%    column 9: Northing in m of LIDAR
%    column 10: Easting in m of LIDAR
%    column 11: heading in rad of LIDAR (North defines 0 rad)
%    column 12: Velocity in m/s of LIDAR
% - map:        image of map in matrix format
% - UTM_map_x:  grid vector in x (column) direction according to xgv in
% function meshgrid
% - UTM_map_y:  grid vector in y (row) direction according to ygv in
% function meshgrid
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
% $Revision: 2.4$ $Date: 2013/05/16$
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Input check
if nargin < 5 || nargin == 6 || nargin > 7
    error('Wrong number of input arguments.')
end
if ~ismatrix(test_series)
    error('test_series must be matrix.')
end
if ~isvector(UTM_map_x) || ~isvector(UTM_map_y)
    error('UTM_map_x and UTM_map_y must be vectors.')
end

if ~exist('usegradient','var')
    usegradient = false;
end

% Global variables
color_code_file=load(color_code);    % Reference file for colour coding
opacity=0.99;                        % Opacity of final image
SNR_limit=1.01;                      % Filter limit of doppler intensity

% Define location to save figures
if exist('root_superposer.txt','file')
    import = importdata('root_superposer.txt');
    root_superposer = cell2mat(import(1,1));
else
    % UI for selecting the data folders
    disp('Please select the root folder for superposes figures and images.')
    root_superposer = uigetdir('','Please select the root folder for figures.');
    
    % Save path to txt-file
    fileID = fopen('root_superposer.txt','wt');
    fprintf(fileID, '%s\n', root_superposer);
    fclose(fileID);
end

% Check existence of folder figures. If not present, create one
root_function = pwd;
latest = 1;
latestmov = 1;
if exist(root_superposer,'dir')~=7
    mkdir(root_superposer);
    disp('Folder "superposed" created"')
elseif size(dir(root_superposer),1) > 2   % Check wether folder is empty
    % Check numbering of figures. If figures are already present, the number
    % should increase by one number
    cd(root_superposer);
    list = dir('*_sequence.fig');
    list = struct2cell(list);
    list = list(1,:);           % Crop list to filenames
    list = char(list);
    A = NaN(size(list,1),1);
    for k = 1:size(list,1)
        temp = strsplit(list(k,:), {'_', '.'},'CollapseDelimiters',true);
        A(k) = str2num(cell2mat(temp(1)));
    end
    latest = max(A) + 1;
    if isempty('latest')
        error('Existing figures could not be identified.')
    end
    % Apply the same to the movie
    list = dir('*_completesequence.avi');
    list = struct2cell(list);
    list = list(1,:);           % Crop list to filenames
    list = char(list);
    A = NaN(size(list,1),1);
    for k = 1:size(list,1)
        temp = strsplit(list(k,:), {'_', '.'},'CollapseDelimiters',true);
        A(k) = str2num(cell2mat(temp(1)));
    end
    latestmov = max(A) + 1;
    if isempty(latest)
        warning('Existing figures could not be identified. Starting with number 1.')
        latest = 1;
    end
    if isempty(latestmov)
        warning('Existing movies could not be identified. Starting with number 1.')
        latestmov = 1;
    end
    cd(root_function)
end

disp(['Figures will be saved beginning with number ', num2str(latest), '.']);
disp(['Movie will be saved with number ', num2str(latestmov), '.']);

% Filter data according to SNR_limit
test_series = test_series((test_series(:,4) >= SNR_limit), :);
% Filter data with range below 50 metres and above 400 metres
test_series = test_series((test_series(:,2) >= 50),:);
% test_series = test_series((test_series(:,2) <= 600),:);

% Filter data for NaN
test_series = test_series(~any(isnan(test_series),2),:);

% Calculate the absolute heading of each beam
head_abs = deg2rad(test_series(:,5) + test_series(:,11));

% Compute the measured times for beam allocation
time = unique(test_series(:,1));

% Calculate the minimal and maximal speed
speed = test_series(:,3)- cos(deg2rad(test_series(:,5))).*test_series(:,12);

if usegradient && exist('limit','var')
    speed_min = 0;
    speed_max = abs(limit);
elseif usegradient && ~exist('limit','var')
    speed_min = min(speed,[],1);
    speed_max = max(speed,[],1);
elseif ~usegradient && exist('limit','var')
    speed_min = -abs(limit);
    speed_max = abs(limit);
elseif ~usegradient && ~exist('limit','var')
    speed_min = min(speed,[],1);
    speed_max = max(speed,[],1);
end

% Create new matrix test_final with reduced entries. The coordinates will
% be changed to x and y coordinates and added to the existing UTM data
%    Column 1:  Time in seconds, starting at 12AM of the measurement day
%    column 2:  absolute radial velocity in m/s (variable speed)
%    column 3:  azimuthal angle in degrees
%    Column 4:  total northing in m of measurement point
%    Column 5:  total easting in m of measurement point
%    Column 6:  Northing in m of LIDAR
%    Column 7:  Easting in m of LIDAR
%    Column 8:  absolute heading of beams in rad
%    Column 9:  absolute speed of a measurement point.

easting_tot = test_series(:,10) + test_series(:,2).*sin(head_abs)...
    .*cos(deg2rad(test_series(:,6)));
northing_tot = test_series(:,9) + test_series(:,2).*cos(head_abs)...
    .*cos(deg2rad(test_series(:,6)));
test_final = [test_series(:,[1 3 5]),...
    northing_tot,easting_tot,...
    test_series(:,[9 10]), head_abs, speed];

% Set up grid with the vector inputs UTM_map_x, UTM_map_y
[x_mesh, y_mesh] = meshgrid(UTM_map_x, UTM_map_y);

% Empty array for plot savings
% v_last = [];

for k=2:length(time)
    test_case1 = test_final(time(k-1) == test_final(:,1),:);
    test_case2 = test_final(time(k) == test_final(:,1),:);
    % Sum the values of the 2 beams and interpolate between them
    x_points = [test_case1(:,5); test_case2(:,5)];
    y_points = [test_case1(:,4); test_case2(:,4)];
    % Check if data is within the range
    if min(x_points,[],1) < x_mesh(1,1) || max(x_points,[],1) > x_mesh(end,end)
        warning('The points in x are outside the range of the map.')
    end
    if min(y_points,[],1) < y_mesh(1,1) || max(y_points,[],1) > y_mesh(end,end)
        warning('The points in y are outside the range of the map.')
    end
    
    speed_data = [test_case1(:,9); test_case2(:,9)];
    
    % Create matrix with dataset of properties to be interpolated (northing,
    % easting, heading, velocity)
    % ts_interp = test_series(isnan(test_series(:,2)),:);
    vr_interp = unique([x_points, y_points, speed_data],'rows');
    
    % For interplation, the function unique is not sufficient to delete enough
    % entries. The time data will be checked for differences of 0, which
    % depends on the machine precision.
    % For x-direction
    delta = diff(vr_interp(:,1));
    nodiff = find(delta == 0);
    for n=1:length(nodiff)
        vr_interp(nodiff(n),1) = NaN;
    end
    vr_interp = vr_interp(~isnan(vr_interp(:,1)),:);
    % For y-direction
    delta = diff(vr_interp(:,2));
    nodiff = find(delta == 0);
    for n=1:length(nodiff)
        vr_interp(nodiff(n),2) = NaN;
    end
    vr_interp = vr_interp(~isnan(vr_interp(:,2)),:);
    
    % Reduce mesh to interpolation area
    [~, x_min] = find(x_mesh(1,:) >= min(vr_interp(:,1)),1,'first');
    [y_min, ~] = find(y_mesh(:,1) >= min(vr_interp(:,2)),1,'first');
    [~, x_max] = find(x_mesh(1,:) <= max(vr_interp(:,1)),1,'last');
    [y_max, ~] = find(y_mesh(:,1) <= max(vr_interp(:,2)),1,'last');
    x_min = x_min-1;
    y_min = y_min-1;
    x_max = x_max+1;
    y_max = y_max+1;
    
    % Check the values for x_min, y_min, x_max, y_max
    if x_min < 1 || isempty(x_min)
        x_min = 1;
    end
    if y_min < 1 || isempty(y_min)
        y_min = 1;
    end
    if isempty(x_max) || x_max > length(UTM_map_x)
        x_max = length(UTM_map_x);
    end
    if isempty(y_max) || y_max > length(UTM_map_y)
        y_max = length(UTM_map_y);
    end
    
    x_mesh_red = x_mesh(y_min:y_max,x_min:x_max);
    y_mesh_red = y_mesh(y_min:y_max,x_min:x_max);
    
    % Compute the interpolation in the reduced area
    v_r = griddata(vr_interp(:,1),vr_interp(:,2),vr_interp(:,3),...
        x_mesh_red, y_mesh_red,'natural');
    
    % Implement data into map conform format
    v_final = NaN(size(map,1),size(map,2));
    
    if size(v_final(y_min:y_max,x_min:x_max)) == size(v_r)
    v_final(y_min:y_max,x_min:x_max) = v_r;
    end
    
    % Compute gradient to enlarge contrast in final plot
    [v_final_x, v_final_y] = gradient(v_final);
    grad_final = sqrt(v_final_x.^2 + v_final_y.^2);
    
    % TESTWRITE FOR USE OF GRADIENT
    if usegradient
    v_final = grad_final;
    end
    
    % Compute the colouring of the map according to radial speed
    color_map = overlay(map, opacity, color_code_file,speed_min,...
        speed_max, v_final);
    
    % Implementing former interpolations into the latest plot
    if k<=40
        v_last(:,:,k-1) = v_final;
        opacity_temp = opacity;
        if k>2
            for kk=k-2:-1:1
                opacity_temp = (1-0.025)*opacity_temp;
                color_map = overlay(color_map, opacity_temp,...
                    color_code_file,speed_min, speed_max, v_last(:,:,kk));
            end
        end
    else
        v_last(:,:,1:38) = v_last(:,:,2:39);
        v_last(:,:,39) = v_final;
        opacity_temp = opacity;
        for kk=39:-1:1
            opacity_temp = (1-0.025)*opacity_temp;
            color_map = overlay(color_map, opacity_temp,...
                color_code_file,speed_min, speed_max, v_last(:,:,kk));
        end
    end

    % Plot the figure
    close all
    h = figure('visible','off');
    imshow(color_map);
    text(30,30,...
        ['time in s: ',num2str(test_case2(1,1))],...
        'Fontsize',16,...
        'BackgroundColor',[1 1 1]);
    colormap jet
    colorbar('YTickLabel',speed_min:(speed_max-speed_min)/6:...
        speed_max,'Fontsize',14);
    if usegradient
        title('magnitude of wind speed gradient (shear) in 1/s','Fontsize',18)
    else
        title('Radial wind speed in m/s','Fontsize',18)
    end
    cd(root_superposer)
    saveas(gcf,sprintf('%04d%s',(k-2+latest), '_sequence.fig'));
    set(gcf,'PaperPositionMode','auto')
    if ~exist('./images','dir')
        mkdir('./images')
    end
    cd('images')
    print('-dpng',...   % print png file
        sprintf('%04d%s',(k-2+latest), '_sequence.png'))
    % mov(k-1) = im2frame(zbuffer_cdata(gcf));
    cd ..
    cd(root_function)
    close all
end
close all

% % Create a movie
% cd(root_superposer)
% vidObj = VideoWriter(sprintf('%03d%s',(latestmov),...
%     '_completesequence.avi'));
% vidObj.FrameRate=4;
% open(vidObj);
% writeVideo(vidObj,mov);
% cd(root_function)

end

function cdata = zbuffer_cdata(hfig)
% Get CDATA from hardcopy using zbuffer

% Need to have PaperPositionMode be auto
orig_mode = get(hfig, 'PaperPositionMode');
set(hfig, 'PaperPositionMode', 'auto');

cdata = hardcopy(hfig, '-Dzbuffer', '-r0');

% Restore figure to original state
set(hfig, 'PaperPositionMode', orig_mode);
end

function color_map = overlay(map, opacity, color_code_file,speed_min,...
    speed_max, v_final)
color_map = map;
for l=1:size(map,1)
    for m=1:size(map,2)
        if isnan(v_final(end+1-l,m))
            color_map(l,m,1) = map(l,m, 1);
            color_map(l,m,2) = map(l,m, 2);
            color_map(l,m,3) = map(l,m, 3);
        else
            pos = ceil(64*(v_final(end+1-l,m)-speed_min)/...
                (speed_max-speed_min));
            if pos <= 0
                pos = 1;
            end
            if pos > size(color_code_file,1)
               pos = size(color_code_file,1);
            end
            color_vec=color_code_file(pos,:);
            color_map(l,m,1) = (1-opacity)*map(l,m, 1) + opacity*uint8(255*color_vec(1));
            color_map(l,m,2) = (1-opacity)*map(l,m, 2) + opacity*uint8(255*color_vec(2));
            color_map(l,m,3) = (1-opacity)*map(l,m, 3) + opacity*uint8(255*color_vec(3));
        end
    end
end
end
