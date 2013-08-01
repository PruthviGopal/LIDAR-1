function [map, UTM_map_x, UTM_map_y, map2real] = map_calibration_v2(map_file,lat,long,scale)
%% function map_calibration_v2
% function [UTM_map_x, UTM_map_y, map2real] = map_calibration_v2(map_file,long,lat,scale)
%
% DESCRIPTION
% The function takes the location of a given map with a reference point
% (lat, long) and computes the coordinates. It returns the imported map and
% 2 grid vectors which are in UTM coordinates.
%
% INPUT
% - map_file:   location of image map
% - lat:        latitude of reference point on map
% - lon:        longitude of reference point of map
% - scale:      scale as shown on map_file
%
% OUTPUT
% - map:        image of map in matrix format
% - UTM_map_x:  grid vector in x (column) direction according to xgv in
% function meshgrid
% - UTM_map_y:  grid vector in y (row) direction according to ygv in
% function meshgrid
% map2real:     Relation between real m and pixel on map
%
% Code by: Mohsen Zendehbad
% Edited by: Markus Schmidt
%
% $Revision: 2.2$ $Date: 2013/04/23 $
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Input check
if nargin ~= 4
    error('Wrong number of input arguments.')
end

if ~ischar(map_file)
    error('map_file must be string.')
end

% Read in the image. Ask the user to define 3 points: 1) the reference
% point according to lat and lon of input. 2) One end of the scale 3) The
% 2nd end of the scale. 2) and 3) must correspond to the input variable
% scale.
map=imread(map_file);
imshow(map);
disp('CLICK ON REFERENCE POINT; THEN ON TWO ENDS OF SCALE LINE');
[x,y]=ginput(3);
[y0, x0] = GPS2UTM(lat,long);

% Calculate the scale of the image. The result map2real is in m/pixel
map2real=scale/sqrt( (x(2)-x(3))^2 + (y(2)-y(3))^2 );

% % Export as a txt file
% to_be_saved(1,:)=[long lat];
% to_be_saved(2,:)=UTM_ref;
% to_be_saved(3,:)=[y(1) x(1)];
% to_be_saved(4,:)=[map2real 0];
% save([map_file(1:end-4) '.txt'],'to_be_saved','-ascii')

% Set up a grid with the total UTM coordinates for the map. The y-vector
% has to be flipped due to the different orientation of map(up-down) and
% UTM(down-up).
UTM_map_x = ((1:size(map,2)) - x(1) )*map2real + x0;
UTM_map_y = ((1:size(map,1)) - (size(map,1)-y(1)+1) )*map2real + y0;

close all;
