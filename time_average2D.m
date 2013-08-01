function [ vel_comps ] = time_average2D(gridprops, range, save_figure,...
    lidar_sync,ldm_sync, vlimit)
%% function time_average2D
% [ vel_comps ] = time_average2D(gridprops, range)
% [ vel_comps ] = time_average2D(gridprops, range, save_figure)
% [ vel_comps ] = time_average2D(gridprops, range, save_figure,...
%     lidar_sync,ldm_sync)
% [ vel_comps ] = time_average2D(gridprops, range, save_figure,...
%     lidar_sync,ldm_sync, vlimit)
% 
% DESCRIPTION
% The function applies a time averaging to the data selected during the
% algorithm. Based on the imported data, the user has to choose to time
% intervals which will be declared as measurement one and two. As output it
% provides a struct with velocity components dependend on their
% interpolation method. In addition, plots will be saved to subfolder
% /figures if save_figure is TRUE.
%
% INPUT
% - gridprops: array with properties of final grid: [xmin, xmax, ymin, ymax, stepsize] 
% - range: radial distance between 2 measurements. Is used to resolve the
% range gate of the LIDAR system
% - save_figure: boolean. If TRUE, figures will be saved. Default is FALSE
% - lidar_sync: sync time in s to add to the data in LIDAR and location
% - ldm_sync: sync time in s to add to the data in LDM
% - vlimit: absolute value of maximal wind speed to be shown. Other values
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
% $Revision: 2.0$ $Date: 2013/05/13 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Global variables
known_times = true;     % Times of measurements are known

% Input Check
if nargin > 6 || nargin < 2
    error('Incorrect number of input arguments')
end

if numel(gridprops) ~= 5
    error('Dimension of gridprops is wrong.')
end

if ~exist('save_figure','var')
    save_figure = false;
end

if ~exist('lidar_sync','var')
    lidar_sync = 0;
end

if ~exist('ldm_sync','var')
    ldm_sync = 0;
end

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

%% Computation of 2D velocity field 
if exist('vlimit','var')
    vel_comps = velocity2D( measure_1, measure_2, gridprops, save_figure, vlimit );
else
    vel_comps = velocity2D( measure_1, measure_2, gridprops, save_figure);
end
end

