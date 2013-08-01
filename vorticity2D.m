function [ omega ] = vorticity2D( varargin )
%% function vorticity2D
% a) function [ omega ] = vorticity2D()
% b) function [ omega ] = vorticity2D(x_mesh, y_mesh, u, v, save_figure, omega_limit )
% c) function [ omega ] = vorticity2D( vel_comps, save_figure, omega_limit )
% 
% DESCRIPTION
% The function has 3 possible input combinations. Option a) is without an
% input argument. An assistant is started to find .txt-files with suitable
% data according to the requirements of 'txt_files' from option b).
%
% Option b) provides an char array with the locations of the .txt-files. In
% addition, parametres for save_figure and omega_limit can be provided.
%
% The option c) takes an existing 2D velocity field, which has an uniform
% grid and interpolated values (with linear, nearest neighbour and natural
% neighbour) and computes the vorticity.
% The numerical derivation is solved with the finite difference method,
% first with a 2nd order centered (2c) and second with a 4th order centered (4c)
% algorithm. Figures can be saved with setting 'save_figure' to true and a
% filtering of the final data via omega_limit.
%
% INPUT b)
% - x_mesh:         x-values in grid format
% - y_mesh:         y-values in grid format
% - u:              horizontal velocity according to grid
% - v:              vertical velocity according to grid
% - save_figure:    boolean. If TRUE, figures will be saved. Default is FALSE
% - omega_limit:    absolute value of maximal vorticity to be shown. Other values
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
% - omega_limit:    absolute value of maximal vorticity to be shown. Other values
% will be cropped and set by NaN in the final plots. If not set, no
% filtering will be applied.
%
% OUTPUT
%   for a) and b) only
%   omega:  vorticity field, 2c, according to input grid
%
%   for c) only
% - omega: struct with data of vorticity
%   x_meshgrid:     x-coordinates in meshgrid format
%   y_meshgrid:     y-coordinates in meshgrid format
%   omega_4c_lin:   vorticity field, 4c, linear interpolation
%   omega_4c_near:  vorticity field, 4c, nearest neighbour
%   omega_4c_nat:   vorticity field, 4c, natural neighbour
%
% Code by: Markus Schmidt
%
% $Revision: 1.1$ $Date: 2013/05/15 $
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
    
    prompt = 'Set limit for omega [none]: ';
    inomega = input(prompt);
    if ~isempty(inomega) && isnumeric(inomega)
        omega_limit = inomega;
    elseif ~isnumeric(inomega)
        error('Input of omega_limit must be numeric.')
    end
    clear infigure inomega
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
% Check wether input is struct
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
    omega_limit = varargin{3};
    end
% Check wether input is char array.
elseif ismatrix(varargin{1})
    x_mesh = varargin{1};
    y_mesh = varargin{2};
    u      = varargin{3};
    v      = varargin{4};
    isimport = true;        % Boolean to decide how many calculations
    if nargin >= 5
        save_figure = varargin{5};
    end
    if nargin == 6
        omega_limit = varargin{6};
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
      omega = derivation_2c(x_mesh,y_mesh,u,v);  
else
    omega_2c_lin = derivation_2c(vel_comps.x_meshgrid, vel_comps.y_meshgrid,...
        vel_comps.u_lin ,vel_comps.v_lin);
    omega_2c_near = derivation_2c(vel_comps.x_meshgrid, vel_comps.y_meshgrid,...
        vel_comps.u_near ,vel_comps.v_near);
    omega_2c_nat = derivation_2c(vel_comps.x_meshgrid, vel_comps.y_meshgrid,...
        vel_comps.u_nat ,vel_comps.v_nat);
end

%% Numerical Derivation with Finite-difference 4th order, central
if ~isimport
    omega_4c_lin = derivation_4c(vel_comps.x_meshgrid, vel_comps.y_meshgrid,...
        vel_comps.u_lin ,vel_comps.v_lin);
    omega_4c_near = derivation_4c(vel_comps.x_meshgrid, vel_comps.y_meshgrid,...
        vel_comps.u_near ,vel_comps.v_near);
    omega_4c_nat = derivation_4c(vel_comps.x_meshgrid, vel_comps.y_meshgrid,...
        vel_comps.u_nat ,vel_comps.v_nat);
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
        list = dir('vorticity2d*');
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

    if isimport
        ssource = {omega};
        sfilename = {'vorticity2d_2c'};
        stitle = {'vorticity in 2D plane (2c)'};
    else
        ssource = {omega_2c_lin, omega_2c_near, omega_2c_nat,...
            omega_4c_lin, omega_4c_near, omega_4c_nat};
        sfilename = {'vorticity2d_lin_2c','vorticity2d_near_2c',...
            'vorticity2d_nat_2c','vorticity2d_lin_4c','vorticity2d_near_4c',...
            'vorticity2d_nat_4c'};
        stitle = {'vorticity in 2D plane (Linear Interpolation, 2c)',...
            'vorticity in 2D plane (Nearest Neighbour, 2c)',...
            'vorticity in 2D plane (Natural Neighbour, 2c)',...
            'vorticity in 2D plane (Linear Interpolation, 4c)',...
            'vorticity in 2D plane (Nearest Neighbour, 4c)',...
            'vorticity in 2D plane (Natural Neighbour, 4c)'};
        x_mesh = vel_comps.x_meshgrid;
        y_mesh = vel_comps.y_meshgrid;
    end
    
    % Save plot entries to struct
    plots = struct('source', ssource, 'filename', sfilename, 'title', stitle);
    
    % Plot the struct
    for k=1:size(plots,2)
        set(gcf,'PaperUnits', 'centimeters', 'PaperType', 'A4',...
            'PaperOrientation','portrait','PaperPositionMode', 'manual',...
            'PaperPosition', [2 1 27.7 21])
        contourf(x_mesh, y_mesh, plots(k).source)
        title(plots(k).title)
        xlabel('horizontal distance to LIDAR in m')
        ylabel('vertical distance to LIDAR in m')
        colorbar('location','SouthOutside')
        colormap('jet')
        if exist('omega_limit','var')
            caxis([-abs(omega_limit) abs(omega_limit)])
        end
        printeps(plots(k).filename);
        savefig(plots(k).filename);
        hfigure(gcf) = figure(gcf+1);
    end   
end
close all

%% Create output struct
% - omega: struct with data of vorticity
%   x_meshgrid:     x-coordinates in meshgrid format
%   y_meshgrid:     y-coordinates in meshgrid format
%   omega_2c_lin:   vorticity field, 2c, linear interpolation
%   omega_2c_near:  vorticity field, 2c, nearest neighbour
%   omega_2c_nat:   vorticity field, 2c, natural neighbour
%   omega_4c_lin:   vorticity field, 4c, linear interpolation
%   omega_4c_near:  vorticity field, 4c, nearest neighbour
%   omega_4c_nat:   vorticity field, 4c, natural neighbour
if ~isimport
    omega = struct(...
        'x_meshgrid', vel_comps.x_meshgrid, 'y_meshgrid', vel_comps.y_meshgrid,...
        'omega_2c_lin', omega_2c_lin, 'omega_2c_near', omega_2c_near,...
        'omega_2c_nat', omega_2c_nat, 'omega_4c_lin', omega_4c_lin,...
        'omega_4c_near', omega_4c_near, 'omega_4c_nat',omega_4c_nat);
else
    % Transform result back to original mesh
    if exist('transposed','var')
        omega = transpose(omega);
    end
    if exist('flipped_ud','var')
        omega = flipud(omega);
    end   
    if exist('flipped_lr','var')
        omega = fliplr(omega);
    end 
end

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
