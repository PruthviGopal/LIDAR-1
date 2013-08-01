function [ measure_1, measure_2 ] = ts_selector( test_series, start1, end1,...
    start2, end2)
%% function ts_selector
% [ measure_1, measure_2 ] = ts_selector( test_series )
% [ measure_1, measure_2 ] = ts_selector( test_series, start1, end1,...
%   start2, end2)
%
% DESCRIPTION
% The function takes an existing imeasurement campaign and extract two
% measurements timewise. The start and end times can be provided as an
% input. If no times are given as an input, the function calls the data as
% a plot and asks the user for selection of the times.
%
% INPUT
% - test_series: matrix with the following data
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
% - start1, end1, start2, end2: numerical arrays of the following format
%               [hh, minmin, ss]
%
% OUTPUT
% - measure_*: matrix with the following data
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
%
% Code by: Markus Schmidt
%
% $Revision: 0.2$ $Date: 2013/05/07 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Input check
if ~(nargin == 1 || nargin == 5)
    error('Wrong number of input arguments.')
end

if nargin == 5
    start1 = 3600*start1(1) + 60*start1(2) + start1(3);
    end1 = 3600*end1(1) + 60*end1(2) + end1(3);
    start2 = 3600*start2(1) + 60*start2(2) + start2(3);
    end2 = 3600*end2(1) + 60*end2(2) + end2(3);
end

%% Show a plot of latitude, longite againgst time to evaluate, which time
% intervals are supposed to be extracted
% Reduce plot to entries with not NaN entries
if nargin == 1
    
    % Plot location to show
    figure(1);
    subplot(1,2,1)
    scatter(test_series(~isnan(test_series(:,9)),1),...
        test_series(~isnan(test_series(:,9)),9))
    xlabel('daytime in seconds')
    ylabel('Northing in m')
    
    subplot(1,2,2)
    scatter(test_series(~isnan(test_series(:,10)),1),...
        test_series(~isnan(test_series(:,10)),10))
    xlabel('daytime in seconds')
    ylabel('Eastin in m')
    
    % Entering the start of measurement 1
    disp('Please give values for times of two measurement intervals. Assume that the wake is between measurement 1 and measurement 2!')
    
    correct = false;
    while ~correct
        prompt_1 = 'Please enter the START of measurement 1(in seconds): ';
        start1 = input(prompt_1);
        
        if start1<test_series(1,1) || start1>test_series(end,1)
            disp('The entered value is outside the measurement range. Please try again!')
        else
            correct = true;
        end
    end
    
    % Entering the END of measurement 1
    correct = false;
    while ~correct
        prompt_2 = 'Please enter the END of measurement 1(in seconds): ';
        end1 = input(prompt_2);
        
        if end1<start1 || end1>test_series(end,1)
            disp('The entered value is outside the allowed range. Please try again!')
        else
            correct = true;
        end
    end
    
    % Entering the START of measurement 2
    correct = false;
    while ~correct
        prompt_3 = 'Please enter the START of measurement 2(in seconds): ';
        start2 = input(prompt_3);
        
        if start2>test_series(end,1)
            disp('The entered value is outside the allowed range. Please try again!')
        else
            correct = true;
        end
    end
    
    % Entering the END of measurement 2
    correct = false;
    while ~correct
        prompt_4 = 'Please enter the END of measurement 2(in seconds): ';
        end2 = input(prompt_4);
        
        if end2<start2 || start2>test_series(end,1)
            disp('The entered value is outside the allowed range. Please try again!')
        else
            correct = true;
        end
    end
    
    % Cleaning the workspace
    clear prompt*
end

%% Crop test_series into defined measurements. If there is no entry with the
% time chosen, the algorithm round up the start and round down the end,
% respectively.

start1_row = find(test_series(:,1) <= start1);
start1_row = start1_row(end);   % round up

end1_row = find(test_series(:,1) <= end1);
end1_row = end1_row(end);

start2_row = find(test_series(:,1) <= start2);
start2_row = start2_row(end);

end2_row = find(test_series(:,1) <= end2);
end2_row = end2_row(end);

measure_1 = test_series(start1_row:end1_row,:);     % Final crop of measurement 1
measure_2 = test_series(start2_row:end2_row,:);     % Final crop of measurement 2

%% Ask the user how WindRover was aligned during measurement
mm = true;
while mm
    prompt_1 = 'Is elevation in both measurements aligned to same direction? [Y/N]: ';
    flip1 = input(prompt_1,'s');
    if flip1 == 'Y' || flip1 =='y'
        disp('Measurement 2 will not be changed.')
        mm = false;
    elseif flip1 == 'N' || flip1 =='n'
        disp('Elevation in measurement 2 will be changed.')
        % Flip the orientation of the measurement
        for n=1:size(measure_2,1)
            measure_2(n,6) = 90 - measure_2(n,6) + 90;
        end
        mm = false;
    else
        disp('Wrong input, please try again!')
    end
end

end

