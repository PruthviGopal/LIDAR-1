function [ phase_vorticity ] = phase_vorticity2D( phase_velocity, save_figure, limit )
%% function phase_vorticity2D
% [ phase_vorticity ] = phase_vorticity2D( phase_velocity )
% [ phase_vorticity ] = phase_vorticity2D( phase_velocity, save_figure )
% [ phase_vorticity ] = phase_vorticity2D( phase_velocity, save_figure, limit )
% 
% DESCRIPTION
% The function takes an input struct phase_velocity, computed by function
% phase_average2D, and calculates the vorticity with a Finite difference
% method of 4th order. If save_figure is set, the function saves figures of
% each result in a subfolder of rootfigures. The input 'limit' can limit
% the shown vorticity range in those figures.
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
% - save_figure: boolean. If TRUE, figures will be saved. Default is FALSE
% - limit: absolute value of maximal vorticity to be shown in plots. If not
% set, no limit will be applied to the final plots.
%
% OUTPUT
% - phase_vorticity: struct with data of the velocity components
%       phase_range: range of phases in rad
%       x_meshgrid: x-coordinates in meshgrid format
%       y_meshgrid: y-coordinates in meshgrid format
%       vorticity_2D_lin: absolute value of wind component (linear)
%       vorticity_2D_near: absolute value of wind component (nearest neighbour)
%       vorticity_2D_nat: absolute value of wind component (natural neighbour)
%
% Code by: Markus Schmidt
%
% $Revision: 0.2$ $Date: 2013/05/15 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

close all

% Global variables
interp_method = 3;          % 1: Linear, 2: Nearest N.,3: Natural N.

% Input Check
if nargin > 3 || nargin == 0
    error('Incorrect number of input arguments')
end

if ~exist('save_figure','var')
    save_figure = false;
end

scriptfolder = pwd;

%% Numerical derivation of vorticity
% Compute number of discrete phases
Nphase = size(phase_velocity,2);

for k=1:Nphase
    % % Numerical Derivation with Finite-difference 2nd order, central
    % omega_2c_lin = derivation_2c(phase_velocity.x_meshgrid, phase_velocity.y_meshgrid,...
    %     phase_velocity.u_lin ,phase_velocity.v_lin);
    % omega_2c_near = derivation_2c(phase_velocity.x_meshgrid, phase_velocity.y_meshgrid,...
    %     phase_velocity.u_near ,phase_velocity.v_near);
    % omega_2c_nat = derivation_2c(phase_velocity.x_meshgrid, phase_velocity.y_meshgrid,...
    %     phase_velocity.u_nat ,phase_velocity.v_nat);
    
    % Numerical Derivation with Finite-difference 4th order, central
    omega_4c_lin = derivation_4c(phase_velocity(k).x_meshgrid, phase_velocity(k).y_meshgrid,...
        phase_velocity(k).u_lin ,phase_velocity(k).v_lin);
    omega_4c_near = derivation_4c(phase_velocity(k).x_meshgrid, phase_velocity(k).y_meshgrid,...
        phase_velocity(k).u_near ,phase_velocity(k).v_near);
    omega_4c_nat = derivation_4c(phase_velocity(k).x_meshgrid, phase_velocity(k).y_meshgrid,...
        phase_velocity(k).u_nat ,phase_velocity(k).v_nat);
    
    % Save results
    phase_vorticity(1,k) = struct('phase_range', phase_velocity(k).phase_range,...
         'x_meshgrid', phase_velocity(k).x_meshgrid,...
         'y_meshgrid', phase_velocity(k).y_meshgrid,...
         'vorticity_2D_lin',omega_4c_lin, 'vorticity_2D_near', omega_4c_near,...
         'vorticity_2D_nat',omega_4c_nat);
end

%% Visualisation of results
% Check existence of folder figures. If not present, create one
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
    latest = 1;
    if exist(rootfigures,'dir')~=7
        mkdir(rootfigures);
    end  
    cd(rootfigures);
    list = dir('phase_vorticity-*');
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
    currentfolder = sprintf('%s-%02d','phase_vorticity',latest);
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
            titlestart = 'vorticity in s^{-1} in 2D plane (Linear Interpolation) for ';
            namestart = 'phase_vort_lin';
            vorticity = NaN(size(phase_vorticity(1).x_meshgrid,1),...
                size(phase_vorticity(1).x_meshgrid,2),Nphase);
            for m=1:Nphase
                if ~isempty(phase_vorticity(m).vorticity_2D_lin)
                    vorticity(:,:,m) = phase_vorticity(m).vorticity_2D_lin;
                else
                    vorticity(:,:,m) = [];
                end
            end
        case 2
            titlestart = 'vorticity in s^{-1} in 2D plane (Nearest N.) for ';
            namestart = 'phase_vort_near';
            vorticity = NaN(size(phase_vorticity(1).x_meshgrid,1),...
                size(phase_vorticity(1).x_meshgrid,2),Nphase);
            for m=1:Nphase
                if ~isempty(phase_vorticity(m).vorticity_2D_near)
                    vorticity(:,:,m) = phase_vorticity(m).vorticity_2D_near;
                else
                    vorticity(:,:,m) = [];
                end
            end
        case 3
            titlestart = 'vorticity in s^{-1} in 2D plane (Natural N.) for ';
            namestart = 'phase_vort_nat';
            vorticity = NaN(size(phase_vorticity(1).x_meshgrid,1),...
                size(phase_vorticity(1).x_meshgrid,2),Nphase);
            for m=1:Nphase
                if ~isempty(phase_vorticity(m).vorticity_2D_nat)
                    vorticity(:,:,m) = phase_vorticity(m).vorticity_2D_nat;
                else
                    vorticity(:,:,m) = [];
                end
            end
    end
    for l=1:Nphase
        if ~isempty(vorticity(:,:,l))
        hfigure(gcf) = figure(gcf+1);
        contourf(phase_vorticity(l).x_meshgrid,...
            phase_vorticity(l).y_meshgrid,vorticity(:,:,l))
        title([titlestart,...
            num2str(phase_vorticity(l).phase_range(1),3), ' to ',...
            num2str(phase_vorticity(l).phase_range(2),3), 'radians.']);
        xlabel('horizontal distance to LIDAR in m')
        ylabel('vertical distance to LIDAR in m')
        if exist('limit','var')
            caxis([-abs(limit) abs(limit)])
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

function [ omega_2c ] = derivation_2c(x_mesh, y_mesh, u ,v)
% Set up matrix for numerical derivation
delta_x = x_mesh(1,2)- x_mesh(1,1);
delta_y = y_mesh(2,1)- y_mesh(1,1);

% Calculate derivaties du/dy and dv/dx
uy_2c = NaN(size(x_mesh));
vx_2c = NaN(size(x_mesh));

%Fill matrix with values
% Formula for algorithm: u' = (-u_(i-1) + u_(i+1))/(2*delta_X)
for k=2:size(u,1)-1         % Borders cannot be calculated
    for l=1:size(u,2)
        uy_2c(k,l) = (-u(k-1,l) + u(k+1,l))/(2*delta_y);
    end
end

for k=1:size(u,1)
    for l=2:size(u,2)-1     % Borders cannot be calculated
        vx_2c(k,l) = (-v(k,l-1) + v(k,l+1))/(2*delta_x);
    end
end
% Output
omega_2c = -uy_2c + vx_2c;
end

function [ omega_4c ] = derivation_4c(x_mesh, y_mesh, u ,v)
% Set up matrix for numerical derivation
delta_x = x_mesh(1,2)- x_mesh(1,1);
delta_y = y_mesh(2,1)- y_mesh(1,1);

% Calculate derivaties du/dy and dv/dx
uy_4c = NaN(size(x_mesh));
vx_4c = NaN(size(x_mesh));

%Fill matrix with values
% % Formula for algorithm:
% u' = (u_(i-2) - 8*u_(i-1)+8*u_(i+1)-u_(j+2))/(12*delta_X)
for k=3:size(u,1)-2         % Borders cannot be calculated
    for l=1:size(u,2)
        uy_4c(k,l) = (u(k-2,l) - 8*u(k-1,l) + 8*u(k+1,l)- u(k+2,l))/(12*delta_y);
    end
end

for k=1:size(u,1)
    for l=3:size(u,2)-2     % Borders cannot be calculated
        vx_4c(k,l) = (v(k,l-2) - 8*v(k,l-1) + 8*v(k,l+1)- v(k,l+2))/(12*delta_x);
    end
end

% Output
omega_4c = -uy_4c + vx_4c;
end
