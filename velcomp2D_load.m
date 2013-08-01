function [ x_mesh, y_mesh, u, v ] = velcomp2D_load( )
%% function velcomp2D_load
% function [ vel_comps ] = velcomp2D_load( )
% 
% DESCRIPTION
% THe function assumed .txt. file for the coordinates and the velocity
% components in MATLAB's meshgrid format. A struct will be created as
% output for further computation.
%
% INPUT
% - none
%
% OUTPUT
%   x_mesh:     x-coordinates in meshgrid format
%   y_mesh:     y-coordinates in meshgrid format
%   u:          horizontal velocity component
%   v:          vertical velocity component
%
% Code by: Markus Schmidt
%
% $Revision: 0.2$ $Date: 2013/04/09 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Error messages
if nargin ~= 0
    error('Incorrect number of input arguments.')
end

% Check with operating system is used. Important to create full path names.
if isunix || ismac
    slash = '/';
else
    slash = '\';
end

% UI for selecting the data files
disp('Please select the file with x-data.')
[xName,xPath,~] = uigetfile('*.txt','Please select the file with x-data:');
x_file = sprintf('%s%s%s',xPath, slash, xName);

disp('Please select the file with y-data.')
[yName,yPath,~] = uigetfile('*.txt','Please select the file with y-data:');
y_file = sprintf('%s%s%s',yPath, slash, yName);

disp('Please select the file with u-data.')
[uName,uPath,~] = uigetfile('*.txt','Please select the file with u-data:');
u_file = sprintf('%s%s%s',uPath, slash, uName);

disp('Please select the file with v-data.')
[vName,vPath,~] = uigetfile('*.txt','Please select the file with v-data:');
v_file = sprintf('%s%s%s',vPath, slash, vName);

% Load the files into MATLAB
x_mesh = load(x_file);
y_mesh = load(y_file);
u = load(u_file);
v = load(v_file);
end

