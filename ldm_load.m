function ldmmat = ldm_load( ldm_data )
%% function ldm_load
% function ldmmat = ldm_load( ldm_data )
% 
% DESCRIPTION
% The function import the data recorder by LDM. The data structure starts
% at 12AM and records the time between two blades passing the LDM.
% In addition, the function looks for missing measurements. If the time
% values changes by more than defined in variable 'threshold', a linear
% interpolation takes place and inserts a new time value between the two
% others. A command output tells how many values have been inserted.
% 
% INPUT
% - ldm_data: .txt file with information of blade passage from LDM. The
% time intervalls are given in miliseconds
% 
% OUTPUT
% - ldmmat: numerical array with one column. The entries show the blade's
% passing time in seconds, starting at 12AM.
%
% Code by: Markus Schmidt
%
% $Revision: 1.1$ $Date: 2013/03/25 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Input check
if nargin~=1
    error('Wrong number of input arguments.')
end

% Import of file:
ldmmat = dlmread(ldm_data);
% Example-value: 3.0269560e+04

% Identify outlines
% Assuming a steady change in number of revolutions, the variance between
% the average of 2 nor given at t_1 and t_3 should not bigger than a
% defined threshold compared to the computed nor at t2. Otherwise a missing
% measurement is assumed.

threshold = 0.01; % definded relative deviation
n_dev = 0; % Number of values inserted due to missing measurements.

for j=2:length(ldmmat)-1
    % The value 'dev' deviation checks the to time intervals. If the
    % deviation is bigger than defined in with threshold, a new value will
    % be inserted.
    dev = ldmmat(j+1)/ldmmat(j);
    if (dev >= 1 + threshold) || (dev <= 1 - threshold)
        ldmmat_av = (ldmmat(j+1) + ldmmat(j))/2;
        ldmmat = [ldmmat(1:j)', ldmmat_av, ldmmat(j+1:end)']';
        n_dev = n_dev + 1;
        j=2; % Restart for-loop
    end
end

% Display the number of interpolated values:
disp(['Number of interpolated time values for LMD import: ', num2str(n_dev)])

end
