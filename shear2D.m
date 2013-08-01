function [ shear_x, shear_y ] = shear2D( varargin )
%% function shear2D
% a) function [ shear_x, shear_y ] = shear2D()
% b) function [ shear_x, shear_y ] = shear2D(x_mesh, y_mesh, u, v, save_figure, shear_limit )
% c) function [ shear_x, shear_y ] = shear2D( vel_comps, save_figure, shear_limit )
% 
% DESCRIPTION
% The function has 3 possible input combinations. Option a) is without an
% input argument. An assistant is started to find .txt-files with suitable
% data according to the requirements of 'txt_files' from option b).
%
% Option b) provides an char array with the locations of the .txt-files. In
% addition, parametres for save_figure and shear_limit can be provided.
%
% The option c) takes an existing 2D velocity field, which has an uniform
% grid and interpolated values (with linear, nearest neighbour and natural
% neighbour) and computes the wind shear.
% The numerical derivation is solved with the finite difference method,
% with a 4th order centered (4c) algorithm. Figures can be saved with setting 'save_figure' to true and a
% filtering of the final data via shear_limit.
%
% INPUT b)
% - x_mesh:         x-values in grid format
% - y_mesh:         y-values in grid format
% - u:              horizontal velocity according to grid
% - v:              vertical velocity according to grid
% - save_figure:    boolean. If TRUE, figures will be saved. Default is FALSE
% - shear_limit:    absolute value of maximal shear to be shown. Other values
% will be cropped and set by NaN in the final plots. If not set, no
% filtering will be applied.
%
% INPUT c)
% - vel_comps:      struct with data of the velocity components
%   x_meshgrid:     x-coordinates in meshgrid format
%   y_meshgrid:     y-coordinates in meshgrid format
%   u_lin:          horizontal wind component (linear)
%   v_lin:          vertical wind component (linear)
%   u_near:         horizontal wind component (nearest neighbour)
%   v_near:         vertical wind component (nearest neighbour)
%   u_nat:          horizontal wind component (natural neighbour)
%   v_nat:          vertical wind component (natural neighbour)
% - save_figure:    boolean. If TRUE, figures will be saved. Default is FALSE
% - shear_limit:    absolute value of maximal shear to be shown. Other values
% will be cropped and set by NaN in the final plots. If not set, no
% filtering will be applied.
%
% OUTPUT
%   shear_x: wind shear in x-direction
%   shearyx: wind shear in y-direction
%
% Code by: Markus Schmidt
%
% $Revision: 0.2$ $Date: 2013/04/10 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Global settings
close all                               % Close all open figures

%% Input check
if nargin>6
    error('Number of input arguments wrong.')
end
 
% Check wether input is present. If not, start import function
if nargin ==0
    [ x_mesh, y_mesh, u, v ] = velcomp2D_load();
    
    prompt = 'Do you want figures? Y/N [N]: ';
    infigure = input(prompt,'s');
    if isempty(infigure) || infigure == 'N' || infigure == 'n'
        save_figure = false;
    elseif infigure == 'Y' || infigure == 'y'
        save_figure = true;
    else
        error('Input of save_figure must be boolean.')
    end
    
    prompt = 'Set limit for shear [none]: ';
    inshear = input(prompt);
    if ~isempty(inshear) && isnumeric(inshear)
        shear_limit = inshear;
    elseif ~isnumeric(inshear)
        error('Input of shear_limit must be numeric.')
    end
    clear infigure inshear
    isimport = true;        % Boolean to decide how many calculations
    
    % Check if grid is provided according to requirements in meshgrid
    % Check if ngrid is present
    if x_mesh(1,1) > x_mesh(end,1)  
        x_mesh = flipud(x_mesh);
        u = flipud(u);
        v = flipud(v);
        flipped_ud = true;
    end
    if (x_mesh(1,2)- x_mesh(1,1)) == 0   
        x_mesh = transpose(x_mesh);
        u = transpose(u);
        v = transpose(v);
        transposed = true;
    end
    if y_mesh(1,1) > y_mesh(1,end)  
        y_mesh = fliplr(y_mesh);
        u = fliplr(u);
        v = fliplr(v);
        flipped_lr = true;
    end
    if (y_mesh(2,1)- y_mesh(1,1)) == 0
        y_mesh = transpose(y_mesh);
        if ~exist('transposed','var')
        u = transpose(u);
        v = transpose(v);
        transposed = true;
        end
    end
% Check wether input is char array.
elseif isstruct(varargin{1})
    vel_comps = varargin{1};
    isimport = false;
    
    % Check correct size and alignment of meshgrid
    if size(vel_comps.x_meshgrid) ~= size(vel_comps.y_meshgrid)
        error('Input meshgrids are not of same size.')
    end
    % Check if grid is provided according to requirements in meshgrid
    if (vel_comps.x_meshgrid(1,2)- vel_comps.x_meshgrid(1,1)) == 0   
        temp = vel_comps.x_meshgrid;
        vel_comps.x_meshgrid = vel_comps.y_meshgrid;
        vel_comps.y_meshgrid = temp;
        clear temp
    end
    if nargin >= 2
    save_figure = varargin{2}; 
    end
    if nargin == 3
    shear_limit = varargin{3};
    end
    
% Check wether input is char array.
elseif isnumeric(varargin{1})
    x_mesh = varargin{1};
    y_mesh = varargin{2};
    u      = varargin{3};
    v      = varargin{4};
    isimport = true;        % Boolean to decide how many calculations
    if nargin >= 5
        save_figure = varargin{5};
    end
    if nargin == 6
        shear_limit = varargin{6};
    end
    % Check if grid is provided according to requirements in meshgrid
    % Check if ngrid is present
    if x_mesh(1,1) > x_mesh(end,1)  
        x_mesh = flipud(x_mesh);
        u = flipud(u);
        v = flipud(v);
        flipped_ud = true;
    end
    if (x_mesh(1,2)- x_mesh(1,1)) == 0   
        x_mesh = transpose(x_mesh);
        u = transpose(u);
        v = transpose(v);
        transposed = true;
    end
    if y_mesh(1,1) > y_mesh(1,end)  
        y_mesh = fliplr(y_mesh);
        u = fliplr(u);
        v = fliplr(v);
        flipped_lr = true;
    end
    if (y_mesh(2,1)- y_mesh(1,1)) == 0
        y_mesh = transpose(y_mesh);
        if ~exist('transposed','var')
        u = transpose(u);
        v = transpose(v);
        transposed = true;
        end
    end
% If input is not correct, abort function 
else
    error('Input does not have right format.')
end

% Check and set variable 'save_figure'
if ~exist('save_figure','var')
    save_figure = 0;
end

%% Numerical Derivation with Finite-difference 2nd order, central
if isimport
      [shear_x, shear_y] = derivation_2c(x_mesh,y_mesh,u,v);
      % Calculate magnitude
      shear_mag = sqrt(shear_x.^2 + shear_y.^2);
else
    [shear_2c_lin_x, shear_2c_lin_y] = derivation_2c(vel_comps.x_meshgrid,...
        vel_comps.y_meshgrid,vel_comps.u_lin ,vel_comps.v_lin);
    [shear_2c_near_x, shear_2c_near_y] = derivation_2c(vel_comps.x_meshgrid,...
        vel_comps.y_meshgrid, vel_comps.u_near ,vel_comps.v_near);
    [shear_2c_nat_x, shear_2c_nat_y] = derivation_2c(vel_comps.x_meshgrid,...
        vel_comps.y_meshgrid, vel_comps.u_nat ,vel_comps.v_nat);
    
    % Calculate magnitude
    shear_2c_lin_mag = sqrt(shear_2c_lin_x.^2 + shear_2c_lin_y.^2);
    shear_2c_near_mag = sqrt(shear_2c_near_x.^2 + shear_2c_near_y.^2);
    shear_2c_nat_mag = sqrt(shear_2c_nat_x.^2 + shear_2c_nat_y.^2);
end
%% Numerical Derivation with Finite-difference 4th order, central
if ~isimport
    [shear_4c_lin_x, shear_4c_lin_y] = derivation_4c(vel_comps.x_meshgrid,...
        vel_comps.y_meshgrid,vel_comps.u_lin ,vel_comps.v_lin);
    [shear_4c_near_x, shear_4c_near_y] = derivation_4c(vel_comps.x_meshgrid,...
        vel_comps.y_meshgrid, vel_comps.u_near ,vel_comps.v_near);
    [shear_4c_nat_x, shear_4c_nat_y] = derivation_4c(vel_comps.x_meshgrid,...
        vel_comps.y_meshgrid, vel_comps.u_nat ,vel_comps.v_nat);
    
    % Calculate magnitude
    shear_4c_lin_mag = sqrt(shear_4c_lin_x.^2 + shear_4c_lin_y.^2);
    shear_4c_near_mag = sqrt(shear_4c_near_x.^2 + shear_4c_near_y.^2);
    shear_4c_nat_mag = sqrt(shear_4c_nat_x.^2 + shear_4c_nat_y.^2);
end

%% Plot results
if save_figure
    % Set properties for plot loop
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
    scriptfolder = pwd;
    if exist(rootfigures,'dir')~=7
        mkdir(rootfigures);
    elseif size(dir(rootfigures),1) > 2   % Check wether folder is empty
        % Check numbering of figures. If figures are already present, the number
        % should increase by one number
        cd(rootfigures);
        list = dir('shear2d*');
        list = struct2cell(list);
        list = list(1,:);           % Crop list to filenames
        list = char(list);
        A = NaN(size(list,1),1);
        for k = 1:size(list,1)
            temp = strsplit(list(k,:), {'-', '.'},'CollapseDelimiters',true);
            A(k) = str2num(cell2mat(temp(2)));
        end
        latest = max(A) + 1;
        if isempty('latest')
            error('Existing plots in figures could not be identified.')
        end
        cd(scriptfolder)
    end
    disp(['Vorticity figures will be saved with number ', num2str(latest), '.']);
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
    
    % Set properties for plot loop
    if isimport
        quivx_source = {shear_x};
        quivy_source = {shear_y};
        cont_source = {shear_mag};
        sfilename = {'shear2d_2c'};
        stitle = {'shear in 2D plane (2c)'};
    else
        quivx_source = {shear_2c_lin_x, shear_2c_near_x, shear_2c_nat_x,...
            shear_4c_lin_x, shear_4c_near_x,shear_2c_nat_x};
        quivy_source = {shear_2c_lin_y, shear_2c_near_y, shear_2c_nat_y,...
            shear_4c_lin_y, shear_4c_near_y,shear_2c_nat_y};
        cont_source = {shear_2c_lin_mag, shear_2c_near_mag, shear_2c_nat_mag,...
            shear_4c_lin_mag, shear_4c_near_mag, shear_4c_nat_mag};
        sfilename = {'shear2d_lin_2c','shear2d_near_2c',...
            'shear2d_nat_2c','shear2d_lin_4c','shear2d_near_4c',...
            'shear2d_nat_4c'};
        stitle = {'shear in 2D plane (Linear Interpolation, 2c)',...
            'shear in 2D plane (Nearest Neighbour, 2c)',...
            'shear in 2D plane (Natural Neighbour, 2c)',...
            'shear in 2D plane (Linear Interpolation, 4c)',...
            'shear in 2D plane (Nearest Neighbour, 4c)',...
            'shear in 2D plane (Natural Neighbour, 4c)'};
        x_mesh = vel_comps.x_meshgrid;
        y_mesh = vel_comps.y_meshgrid;
    end
    
    % Save plot entries to struct
    plots = struct('cont_source', cont_source,'quivx_source', quivx_source,...
        'quivy_source', quivy_source,'filename', sfilename, 'title', stitle);
    
    % Plot the struct
    for k=1:size(plots,2)
        set(gcf,'PaperUnits', 'centimeters', 'PaperType', 'A4',...
            'PaperOrientation','portrait','PaperPositionMode', 'manual',...
            'PaperPosition', [2 1 27.7 21])
        contourf(x_mesh, y_mesh, plots(k).cont_source)
        quiver(x_mesh, y_mesh, plots(k).quivx_source,plots(k).quivy_source)
        title(plots(k).title)
        xlabel('horizontal distance to LIDAR in m')
        ylabel('vertical distance to LIDAR in m')
        colorbar('location','SouthOutside')
        colormap('jet')
        if exist('shear_limit','var')
            caxis([0 abs(shear_limit)])
        end
        printeps(plots(k).filename);
        savefig(plots(k).filename);
        hfigure(gcf) = figure(gcf+1);
    end
end
close all

%% Create output

if ~isimport
    shear_x = shear_2c_lin_x;
    shear_y = shear_2c_lin_y;
else
    % Transform result back to original mesh
    if exist('transposed','var')
        shear_x = transpose(shear_x);
        shear_y = transpose(shear_y);
    end
    if exist('flipped_ud','var')
        shear_x = flipud(shear_x);
        shear_y = flipud(shear_y);
    end
    if exist('flipped_lr','var')
        shear_x = fliplr(shear_x);
        shear_y = fliplr(shear_y);
    end   
end

end

function [ shear_x, shear_y] = derivation_2c(x_mesh, y_mesh, u ,v)
% Set up matrix for numerical derivation
delta_x = x_mesh(1,2)- x_mesh(1,1);
delta_y = y_mesh(2,1)- y_mesh(1,1);

% Grid check
if delta_x == 0
    error('Wrong x-mesh as input.')
end

if delta_y == 0
    error('Wrong y-mesh as input.')
end

% Calculate derivaties du/dx, du/dy, dv/dx and dv/dy
ux_2c = NaN(size(x_mesh));
uy_2c = NaN(size(x_mesh));
vx_2c = NaN(size(x_mesh));
vy_2c = NaN(size(x_mesh));

%Fill matrix with values
% Formula for algorithm: u' = (-u_(i-1) + u_(i+1))/(2*delta_X)

% For du/dx
for k=1:size(u,1)
    for l=2:size(u,2)-1     % Borders cannot be calculated
        ux_2c(k,l) = (-u(k,l-1) + u(k,l+1))/(2*delta_x);
    end
end

% For du/dy
for k=2:size(u,1)-1         % Borders cannot be calculated
    for l=1:size(u,2)
        uy_2c(k,l) = (-u(k-1,l) + u(k+1,l))/(2*delta_y);
    end
end

% For dv/dx
for k=1:size(u,1)
    for l=2:size(u,2)-1     % Borders cannot be calculated
        vx_2c(k,l) = (-v(k,l-1) + v(k,l+1))/(2*delta_x);
    end
end

% For dv/dy
for k=2:size(u,1)-1         % Borders cannot be calculated
    for l=1:size(u,2)
        vy_2c(k,l) = (-v(k-1,l) + v(k+1,l))/(2*delta_y);
    end
end

% Output
shear_x = ux_2c + vx_2c;
shear_y = uy_2c + vy_2c;
end

function [ shear_x, shear_y ] = derivation_4c(x_mesh, y_mesh, u ,v)
% Set up matrix for numerical derivation
delta_x = x_mesh(1,2)- x_mesh(1,1);
delta_y = y_mesh(2,1)- y_mesh(1,1);

% Grid check
if delta_x == 0
    error('Wrong x-mesh as input.')
end

if delta_y == 0
    error('Wrong y-mesh as input.')
end

% Calculate derivaties du/dx, du/dy, dv/dx and dv/dy
ux_4c = NaN(size(x_mesh));
uy_4c = NaN(size(x_mesh));
vx_4c = NaN(size(x_mesh));
vy_4c = NaN(size(x_mesh));

%Fill matrix with values
% % Formula for algorithm:
% u' = (u_(i-2) - 8*u_(i-1)+8*u_(i+1)-u_(j+2))/(12*delta_X)

% For du/dx
for k=1:size(u,1)
    for l=3:size(u,2)-2     % Borders cannot be calculated
        ux_4c(k,l) = (u(k,l-2) - 8*u(k,l-1) + 8*u(k,l+1)- u(k,l+2))/(12*delta_x);
    end
end

% For du/dy
for k=3:size(u,1)-2         % Borders cannot be calculated
    for l=1:size(u,2)
        uy_4c(k,l) = (u(k-2,l) - 8*u(k-1,l) + 8*u(k+1,l)- u(k+2,l))/(12*delta_y);
    end
end

% For dv/dx
for k=1:size(u,1)
    for l=3:size(u,2)-2     % Borders cannot be calculated
        vx_4c(k,l) = (v(k,l-2) - 8*v(k,l-1) + 8*v(k,l+1)- v(k,l+2))/(12*delta_x);
    end
end

% For dv/dy
for k=3:size(u,1)-2         % Borders cannot be calculated
    for l=1:size(u,2)
        vy_4c(k,l) = (v(k-2,l) - 8*v(k-1,l) + 8*v(k+1,l)- v(k+2,l))/(12*delta_y);
    end
end

% Output
shear_x = ux_4c + vx_4c;
shear_y = uy_4c + vy_4c;
end
