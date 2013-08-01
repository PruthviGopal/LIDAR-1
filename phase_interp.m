function [ phase_comps ] = phase_interp( measure_phase_1, measure_phase_2,...
    gridprops, phase_range, distance, limit)
%% function phase_interp
% function [ phase_comps ] = phase_interp( measure_phase_1, measure_phase_2,...
%   gridprops, phase_range, distance, limit)
% 
% DESCRIPTION
% The function takes 2 measurements and interpolates a velocity field out
% of them. The measurements will be filtered according to the values in
% phase range. The output can be limited by the value limits. As output one
% receive a struct with the computed velocity fields. If one measurement
% does not contain measurements within the phase, the function returns the
% struct only with the phase range.
%
% INPUT
% - measure_phase_1, measure_phase_2: 
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
%    Column 12: Phase of blade in radians
% - gridprops: array with properties of final grid:
%              [xmin, xmax, ymin, ymax, stepsize] 
% - phase_range: 1x2 numerical array which contains the start and end value
% of the phase to be evaluated. Example phase = [phase_min, phase_max]
%   distance: the distance between the 2 measurements in metres
% - limit: absolute value of maximal wind speed to be shown. Other values
% will be cropped and set by NaN in the final plots. Zero by default
%
% OUTPUT
% - phase_comps: struct with data of the velocity components
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
% $Revision: 0.3$ $Date: 2013/05/06$
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Input check
if nargin < 5 || nargin > 6
    error('Wrong number of input arguments.')
end

% Select the relevant datasets according to phase_range and range large 50
% metres
measure_red_1 = measure_phase_1(measure_phase_1(:,12) >= phase_range(1)...
    & measure_phase_1(:,12) <= phase_range(2),:);
measure_red_2 = measure_phase_2(measure_phase_2(:,12) >= phase_range(1)...
    & measure_phase_2(:,12) <= phase_range(2),:);

% Create output of available measurements
disp(['Phase range ', num2str(phase_range(1)),' to ', num2str(phase_range(2))])
disp(['No. of data in measurement 1: ', num2str(size(measure_red_1,1))])
disp(['No. of data in measurement 2: ', num2str(size(measure_red_2,1))])

% Create desired grid for output
xmin = gridprops(1);
xmax = gridprops(2);
ymin = gridprops(3);
ymax = gridprops(4);
stepsize = gridprops(5);

[x_mesh, y_mesh] = meshgrid(xmin:stepsize:xmax, ymin:stepsize:ymax);

% If no data is present, abort function
if isempty(measure_red_1)
    warning(['Measurement 1 does not contain data between ',...
        num2str(phase_range(1)), ' and ', num2str(phase_range(2)), ' radians.'])
end

if isempty(measure_red_2)
    warning(['Measurement 2 does not contain data between ',...
        num2str(phase_range(1)), ' and ', num2str(phase_range(2)), ' radians.'])
end

if isempty(measure_red_1) || isempty(measure_red_2)
    phase_comps = struct('phase_range',phase_range,'x_meshgrid',x_mesh,...
        'y_meshgrid',y_mesh,'u_lin',[], 'v_lin',[],'velocity_2D_lin', [],...
        'u_near',[], 'v_near',[],'velocity_2D_near',[],...
        'u_nat',[], 'v_nat',[],'velocity_2D_nat', []);
    return
end

% Use function interp_radial2D for interpolation
plane1 = interp_radial2D(measure_red_1);
plane2 = interp_radial2D(measure_red_2);

% Create values from scatteredInterpolant class. Distance is applied to
% measurement to to take the shift into account
vr_lin1 = plane1.linear_interp(x_mesh, y_mesh);
vr_lin2 = plane2.linear_interp(x_mesh-distance, y_mesh);

vr_near1 = plane1.nearest_nb(x_mesh, y_mesh);
vr_near2 = plane2.nearest_nb(x_mesh-distance, y_mesh);

vr_nat1 = plane1.natural_nb(x_mesh, y_mesh);
vr_nat2 = plane2.natural_nb(x_mesh-distance, y_mesh);

% Check wether interpolation is possible
no_interp = isempty(vr_lin1) || isempty(vr_lin2) || ...
    isempty(vr_near1) || isempty(vr_near1) || ...
    isempty(vr_nat1) || isempty(vr_nat2);
if no_interp
    warning(['Interpolation not possible between '...
        num2str(phase_range(1)), ' and ', num2str(phase_range(2)),...
        ' radians.'])
    phase_comps = struct('phase_range',phase_range,'x_meshgrid',x_mesh,...
        'y_meshgrid',y_mesh,'u_lin',[], 'v_lin',[],'velocity_2D_lin', [],...
        'u_near',[], 'v_near',[],'velocity_2D_near',[],...
        'u_nat',[], 'v_nat',[],'velocity_2D_nat', []);
    return
end

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

%% Create output struct
phase_comps = struct('phase_range',phase_range,'x_meshgrid',x_mesh,'y_meshgrid',y_mesh,...
    'u_lin',u_lin, 'v_lin',v_lin,'velocity_2D_lin', velocity_2D_lin,...
    'u_near',u_near, 'v_near',v_near,'velocity_2D_near',velocity_2D_near,...
    'u_nat',u_nat, 'v_nat',v_nat,'velocity_2D_nat', velocity_2D_nat);
end

