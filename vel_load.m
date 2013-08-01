function velmat = vel_load(vel_file, range, correct_yes, heading_angle)
%% function vel_load
% function velmat = vel_load(vel_file, range, correct_yes, heading_angle)
% 
% DESCRIPTION
% This function is used to import the outpout files of the LIDAR system in
% .scn format into Matlab.
% The input variable vel_file is necessarily needed. In addition,
% correct_yes and heading_angle can be given as parameters als well
% 
% INPUT
% - vel_file: The input file in .scn format.
% - range: radial distance between 2 measurements. Is used to resolve the
% range gate of the LIDAR system
% - correct_yes: Boolean. If set to 1 (True) the position will be adjusted
% with the measured values of azimuthal and elevation angle as well as
% heading_angle.
% - heading_angle: Input parameter in degree. Defines the heading angle
% between WindRower and rotational axis of wind turbine. (??)
% 
% OUTPUT
% - vel_mat: numerical array of presicion double. The following data is
% stored:
% column 1: range in m. Radial distance from LIDAR
% column 2: Radial velocity in m/s
% column 3: Doppler intensity
% column 4: Time in seconds, starting at 12AM of the measurement day
% column 5: azimuthal angle in degrees
% column 6: elevation angle in degrees
% column 7: pitch in ??
% column 8: roll in ??
%
% Code by: Mohsen Zendehbad
% Edited by: Markus Schmidt 
%
% $Revision: 1.2$ $Date: 2013/04/05 $

% Error message, if number of input parameters is not correct. Must be 1 or
% 3)
if (nargin~=2) && (nargin~=4)
    error('Wrong number of input parameters in function vel_load. Must be 2 or 4.')
end

% Import of file:
% Independend of number of input parametres, the .scn file is imported
fid=fopen(vel_file);    % Open the 'vel_file' for read access and returns integer file identifier 'fid' 

% Scan 'fid until the entry of Pitch Roll. Afterwards the recorded data
% starts.
b=textscan(fid,'%[^E]');
b=textscan(fid,'%[El]');
b=textscan(fid,'%[PitchRoll]');

% Import the data available in the scn. file, the ray time is divided into
% hours, minutes and seconds by %[:]. Precision is double (parametre %f).
velmat=textscan(fid,'%f %f %f %f %[:] %f %[:] %f %f %f %f %f');

%01 range gate -> 1
%02 velocity -> 2
%03 doppler Intensity -> 3
%04 hour ->
%06 minute ->
%08 second -> 4
%09 azimuthal angle -> 5
%10 elevation angle -> 6
%11 pitch -> 7
%12 roll -> 8

% Convert time into seconds and save in one cell. The range is calculated.
velmat={[(velmat{1}+1)*range velmat{2} velmat{3} velmat{4}*3600+velmat{6}*60+velmat{8} velmat{9} velmat{10} velmat{11} velmat{12}]};
velmat=cell2mat(velmat);    % Convert vel_mat from cell array to numeric array

fclose('all');      % Close all open files

% If the number of input parametres is 3, the correction will be applied
if nargin==4
    if correct_yes==1
        for i=1:size(velmat,1)
            [velmat(i,5),velmat(i,6)]=az_el_reverse_convertor(velmat(i,5),velmat(i,6),heading_angle,velmat(i,7),velmat(i,8));
        end
    end
end
end

%% function reverse_convertor
% function [az_abs,el_abs]=az_el_reverse_convertor(az_rel,el_rel,heading,pitch,roll)
% DESCRIPTION
% 

function [az_abs,el_abs]=az_el_reverse_convertor(az_rel,el_rel,heading,pitch,roll)
    [vec_rel(1,:),vec_rel(2,:),vec_rel(3,:)]=sph2cart(deg2rad(az_rel),deg2rad(el_rel),1);
    
    [y_rolled]=axis_rotation([0;1;0],[1;0;0],-roll);
    [vec_rel_pitched]=axis_rotation(vec_rel,y_rolled,-pitch);
    [vec_rel_pitched_rolled]=axis_rotation(vec_rel_pitched,[1;0;0],-roll);
    
    [az_temp,el_abs,rad]=cart2sph(vec_rel_pitched_rolled(1,:),vec_rel_pitched_rolled(2,:),vec_rel_pitched_rolled(3,:));
    el_abs=rad2deg(el_abs);
    az_temp=rad2deg(az_temp);
    az_abs=az_temp+heading;
    if az_abs<0
        az_abs=az_abs+360;
    else if az_abs>360
            az_abs=az_abs-360;
        end
    end
    el_abs=el_abs';
    az_abs=az_abs';
end

%% function axis rotation
% function [rotated_vec]=axis_rotation(vec,axis,angle)
function [rotated_vec]=axis_rotation(vec,axis,angle)
rotated_vec=...
    [cosd(angle)+axis(1)^2*(1-cosd(angle)) axis(1)*axis(2)*(1-cosd(angle))-axis(3)*sind(angle) axis(1)*axis(3)*(1-cosd(angle))+axis(2)*sind(angle);...
     axis(1)*axis(2)*(1-cosd(angle))+axis(3)*sind(angle) cosd(angle)+axis(2)^2*(1-cosd(angle)) axis(2)*axis(3)*(1-cosd(angle))-axis(1)*sind(angle);...
     axis(1)*axis(3)*(1-cosd(angle))-axis(2)*sind(angle) axis(2)*axis(3)*(1-cosd(angle))+axis(1)*sind(angle) cosd(angle)+axis(3)^2*(1-cosd(angle))]...
     *vec;
end