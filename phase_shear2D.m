function [ phase_shear ] = phase_shear2D( phase_velocity, save_figure, shear_limit )
%% function shear2D
% function [ phase_shear ] = shear2D(phase_velocity)
% function [ phase_shear ] = shear2D(phase_velocity, save_figure)
% function [ phase_shear ] = shear2D(phase_velocity, save_figure, shear_limit)
% 
% DESCRIPTION
% The function takes an input struct phase_velocity, computed
% by function phase_average2D, and calculates the wind shear with a finite
% difference method of 4th order and interpolated values (with linear,
% nearest neighbour and natural neighbour) and computes the wind shear.
%
% INPUT
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
% - save_figure:    boolean. If TRUE, figures will be saved. Default is FALSE
% - shear_limit:    absolute value of maximal shear to be shown. Other values
% will be cropped and set by NaN in the final plots. If not set, no
% filtering will be applied.
%
% OUTPUT
% - phase_shear: struct with data of shear
%       phase_range: range of phases in rad
%       x_meshgrid:     x-coordinates in meshgrid format
%       y_meshgrid:     y-coordinates in meshgrid format
%       shear_2D_lin:   shear field, 4c, linear interpolation
%       shear_2D_lin_x: shear field, 4c, linear interpolation, x-derivative
%       shear_2D_lin_y: shear field, 4c, linear interpolation, y-derivative
%       shear_2D_near:  shear field, 4c, nearest neighbour
%       shear_2D_near_x:shear field, 4c, nearest neighbour, x-derivative
%       shear_2D_near_y:shear field, 4c, nearest neighbour, y-derivative
%       shear_2D_nat:   shear field, 4c, natural neighbour
%       shear_2D_nat_x: shear field, 4c, natural neighbour, x-derivative
%       shear_2D_nat_y: shear field, 4c, natural neighbour, y-derivative
%
% Code by: Markus Schmidt
%
% $Revision: 0.1$ $Date: 2013/05/14 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Global settings
close all                               % Close all open figures

% Global variables
interp_method = 1;          % 1: Linear, 2: Nearest N.,3: Natural N.

%% Input check
if nargin > 3
    error('Number of input arguments wrong.')
end

% Check and set variable 'save_figure'
if ~exist('save_figure','var')
    save_figure = false;
end

%% Compute the shear
% Compute number of discrete phases
Nphase = size(phase_velocity,2);

for k=1:Nphase
    % Numerical Derivation with Finite-difference 2nd order, central
    [shear_2c_lin_x, shear_2c_lin_y] = derivation_2c(phase_velocity(1,k).x_meshgrid,...
        phase_velocity(1,k).y_meshgrid,phase_velocity(1,k).u_lin ,phase_velocity(1,k).v_lin);
    [shear_2c_near_x, shear_2c_near_y] = derivation_2c(phase_velocity(1,k).x_meshgrid,...
        phase_velocity(1,k).y_meshgrid, phase_velocity(1,k).u_near ,phase_velocity(1,k).v_near);
    [shear_2c_nat_x, shear_2c_nat_y] = derivation_2c(phase_velocity(1,k).x_meshgrid,...
        phase_velocity(1,k).y_meshgrid, phase_velocity(1,k).u_nat ,phase_velocity(1,k).v_nat);
    
    % Calculate magnitude
    shear_2c_lin_mag = sqrt(shear_2c_lin_x.^2 + shear_2c_lin_y.^2);
    shear_2c_near_mag = sqrt(shear_2c_near_x.^2 + shear_2c_near_y.^2);
    shear_2c_nat_mag = sqrt(shear_2c_nat_x.^2 + shear_2c_nat_y.^2);
    
    % Numerical Derivation with Finite-difference 4th order, central
    [shear_4c_lin_x, shear_4c_lin_y] = derivation_4c(phase_velocity(1,k).x_meshgrid,...
        phase_velocity(1,k).y_meshgrid,phase_velocity(1,k).u_lin ,phase_velocity(1,k).v_lin);
    [shear_4c_near_x, shear_4c_near_y] = derivation_4c(phase_velocity(1,k).x_meshgrid,...
        phase_velocity(1,k).y_meshgrid, phase_velocity(1,k).u_near ,phase_velocity(1,k).v_near);
    [shear_4c_nat_x, shear_4c_nat_y] = derivation_4c(phase_velocity(1,k).x_meshgrid,...
        phase_velocity(1,k).y_meshgrid, phase_velocity(1,k).u_nat ,phase_velocity(1,k).v_nat);
    
    % Calculate magnitude
    shear_4c_lin_mag = sqrt(shear_4c_lin_x.^2 + shear_4c_lin_y.^2);
    shear_4c_near_mag = sqrt(shear_4c_near_x.^2 + shear_4c_near_y.^2);
    shear_4c_nat_mag = sqrt(shear_4c_nat_x.^2 + shear_4c_nat_y.^2);

    % Save results
    phase_shear(1,k) = struct('phase_range', phase_velocity(k).phase_range,...
        'x_meshgrid', phase_velocity(k).x_meshgrid,...
        'y_meshgrid', phase_velocity(k).y_meshgrid,...
        'shear_2D_lin',shear_4c_lin_mag,...
        'shear_2D_lin_x',shear_4c_lin_x,'shear_2D_lin_y',shear_4c_lin_y,...
        'shear_2D_near',shear_4c_near_mag,...
        'shear_2D_near_x',shear_4c_near_x,'shear_2D_near_y',shear_4c_near_y,...
        'shear_2D_nat',shear_4c_nat_mag,...
        'shear_2D_nat_x',shear_4c_nat_x,'shear_2D_nat_y',shear_4c_nat_y);
end
%% Plot results
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
    end  
    cd(rootfigures);
    list = dir('phase_shear-*');
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
    currentfolder = sprintf('%s-%02d','phase_shear',latest);
    mkdir(currentfolder);
    disp(['Figures will be saved in folder number ', num2str(latest), '.']);
    cd(scriptfolder)
    
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
            titlestart = 'Wind Shear in s^{-1} in 2D plane (Linear Interpolation) for ';
            namestart = 'phase_shear_lin';
            shear = NaN(size(phase_shear(1).x_meshgrid,1),...
                size(phase_shear(1).x_meshgrid,2),Nphase);
            shear_x = NaN(size(phase_shear(1).x_meshgrid,1),...
                size(phase_shear(1).x_meshgrid,2),Nphase);
            shear_y = NaN(size(phase_shear(1).x_meshgrid,1),...
                size(phase_shear(1).x_meshgrid,2),Nphase);
            for m=1:Nphase
                if ~isempty(phase_shear(m).shear_2D_lin)
                    shear(:,:,m) = phase_shear(m).shear_2D_lin;
                    shear_x(:,:,m) = phase_shear(m).shear_2D_lin_x;
                    shear_y(:,:,m) = phase_shear(m).shear_2D_lin_y;
                else
                    shear(:,:,m) = [];
                    shear_x(:,:,m) = [];
                    shear_y(:,:,m) = [];
                end
            end
        case 2
            titlestart = 'Wind Shear in s^{-1} in 2D plane (Nearest N.) for ';
            namestart = 'phase_shear_near';
            shear = NaN(size(phase_shear(1).x_meshgrid,1),...
                size(phase_shear(1).x_meshgrid,2),Nphase);
            shear_x = NaN(size(phase_shear(1).x_meshgrid,1),...
                size(phase_shear(1).x_meshgrid,2),Nphase);
            shear_y = NaN(size(phase_shear(1).x_meshgrid,1),...
                size(phase_shear(1).x_meshgrid,2),Nphase);
            for m=1:Nphase
                if ~isempty(phase_shear(m).shear_2D_near)
                    shear(:,:,m) = phase_shear(m).shear_2D_near;
                    shear_x(:,:,m) = phase_shear(m).shear_2D_near_x;
                    shear_y(:,:,m) = phase_shear(m).shear_2D_near_y;
                else
                    shear(:,:,m) = [];
                    shear_x(:,:,m) = [];
                    shear_y(:,:,m) = [];
                end
            end
        case 3
            titlestart = 'Wind Shear in s^{-1} in 2D plane (Natural N.) for ';
            namestart = 'phase_shear_nat';
            shear = NaN(size(phase_shear(1).x_meshgrid,1),...
                size(phase_shear(1).x_meshgrid,2),Nphase);
            shear_x = NaN(size(phase_shear(1).x_meshgrid,1),...
                size(phase_shear(1).x_meshgrid,2),Nphase);
            shear_y = NaN(size(phase_shear(1).x_meshgrid,1),...
                size(phase_shear(1).x_meshgrid,2),Nphase);
            for m=1:Nphase
                if ~isempty(phase_shear(m).shear_2D_nat)
                    shear(:,:,m) = phase_shear(m).shear_2D_nat;
                    shear_x(:,:,m) = phase_shear(m).shear_2D_nat_x;
                    shear_y(:,:,m) = phase_shear(m).shear_2D_nat_y;
                else
                    shear(:,:,m) = [];
                    shear_x(:,:,m) = [];
                    shear_y(:,:,m) = [];
                end
            end
    end
    for l=1:Nphase
        if ~isempty(shear(:,:,l))
        hfigure(gcf) = figure(gcf+1);
        contourf(phase_shear(l).x_meshgrid,...
            phase_shear(l).y_meshgrid,shear(:,:,l))
        title([titlestart,...
            num2str(phase_shear(l).phase_range(1),3), ' to ',...
            num2str(phase_shear(l).phase_range(2),3), 'radians.']);
        xlabel('horizontal distance to LIDAR in m')
        ylabel('vertical distance to LIDAR in m')
        if exist('shear_limit','var')
            caxis([0 shear_limit])
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
