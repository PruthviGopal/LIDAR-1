function [measure_phase_1, measure_phase_2] = synth_phase(N)
%% function synth_phase
% function [measure_phase_1, measure_phase_2] = synth_phase(N)
% 
% DESCRIPTION The function creates synthetic data which can be used for
% testing inside of function phase_average2D. The input variable defines
% how many phases should be created. The synthetised data is based on a
% potential flow, in which the vorticity gamma1 represent the big swirl of
% the turbine and 3 smaller ones at the tips of the blades. Randomized
% noise is applied to the potential flow in the end.
%
% INPUT
% - N: Number of randomized phases
%
% OUTPUT
% - measure_phase_*: 
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
%
% Code by: Markus Schmidt
%
% $Revision: 0.2$ $Date: 2013/05/07$
%
% This code is licensed under a Creative Commons Attribution-ShareAlike
% 3.0 Unported License
% ( http://creativecommons.org/licenses/by-sa/3.0/deed.en_GB )

% Global parametres
gamma1 = 20;        % Vorticity of turbine
gamma2 = 1;         % Vorticity of blade-end wakes
k      = 5e-40;     % Damping factor for tip vorticities
beta = 0.003;        % Stretching factor of damping
radius = 60;        % Radius of turbine blades
x1     = 10;        % x-position of hub
y1     = 100;       % y-position of hub
distance = 50;      % Distance between 2 measurements

% Grid properties
xmin = -100;   xmax = 100; ymin = 0; ymax = 180; stepwidth = 10;

% Create grid
[x_mesh, y_mesh] = meshgrid(xmin:stepwidth:xmax, ymin:stepwidth:ymax);

% Create random phi
phi = 120.*rand(N,1);

% Create velocity field for each phase
u1  = -gamma1/(2*pi).*(y_mesh-y1)./((x_mesh-x1).^2 + (y_mesh-y1).^2);
u1(isnan(u1)) = 0;
um1 = gamma1/(2*pi).*(y_mesh+y1)./((x_mesh-x1).^2 + (y_mesh+y1).^2);
um1(isnan(um1)) = 0;

v1  = -gamma1/(2*pi).*(x_mesh-x1)./((x_mesh-x1).^2 + (y_mesh-y1).^2);
v1(isnan(v1)) = 0;
vm1 = gamma1/(2*pi).*(x_mesh+x1)./((x_mesh-x1).^2 + (y_mesh+y1).^2);
vm1(isnan(vm1)) = 0;

measure_phase_1 = [];
measure_phase_2 = [];

for k=1:N
    x21 = x1 - radius * cos(deg2rad(phi(k)));
    y21 = y1 - radius * sin(deg2rad(phi(k)));
    x22 = x1 - radius * cos(deg2rad(phi(k) +120));
    y22 = y1 - radius * sin(deg2rad(phi(k) +120));
    x23 = x1 - radius * cos(deg2rad(phi(k) -120));
    y23 = y1 - radius * sin(deg2rad(phi(k) -120));
    
    % Calculation of u
    u21 = (-gamma2/(2*pi).*(y_mesh-y21)./((x_mesh-x21).^2+(y_mesh-y21).^2))...
        *50.*exp(-k*sqrt((x_mesh-x21).^2 +(y_mesh-y21).^2).*beta);
    um21 = (gamma2/(2*pi).*(y_mesh+y21)./((x_mesh-x21).^2+(y_mesh+y21).^2))...
        *50.*exp(-k*sqrt((x_mesh-x21).^2 +(y_mesh+y21).^2).*beta);
    u22 = (-gamma2/(2*pi).*(y_mesh-y22)./((x_mesh-x22).^2+(y_mesh-y22).^2))...
        *50.*exp(-k*sqrt((x_mesh-x22).^2 +(y_mesh-y22).^2).*beta);
    um22 = (gamma2/(2*pi).*(y_mesh+y22)./((x_mesh-x22).^2+(y_mesh+y22).^2))...
        *50.*exp(-k*sqrt((x_mesh-x22).^2 +(y_mesh+y22).^2).*beta);
    u23 = (-gamma2/(2*pi).*(y_mesh-y23)./((x_mesh-x23).^2+(y_mesh-y23).^2))...
        *50.*exp(-k*sqrt((x_mesh-x23).^2 +(y_mesh-y23).^2).*beta);
    um23 = (gamma2/(2*pi).*(y_mesh+y23)./((x_mesh-x23).^2+(y_mesh+y23).^2))...
        *50.*exp(-k*sqrt((x_mesh-x23).^2 +(y_mesh+y23).^2).*beta);
    u = u1 + um1 + u21 + um21 + u22 + um22 + u23 + um23;
    
    % Calculation of v
    v21 = (-gamma2/(2*pi).*(x_mesh-x21)./((x_mesh-x21).^2+(y_mesh-y21).^2))...
        *50.*exp(-k*sqrt((x_mesh-x21).^2 +(y_mesh-y21).^2).*beta);
    vm21 = (gamma2/(2*pi).*(x_mesh+x21)./((x_mesh-x21).^2+(y_mesh+y21).^2))...
        *50.*exp(-k*sqrt((x_mesh-x21).^2 +(y_mesh+y21).^2).*beta);
    v22 = (-gamma2/(2*pi).*(x_mesh-x22)./((x_mesh-x22).^2+(y_mesh-y22).^2))...
        *50.*exp(-k*sqrt((x_mesh-x22).^2 +(y_mesh-y22).^2).*beta);
    vm22 = (gamma2/(2*pi).*(x_mesh+x22)./((x_mesh-x22).^2+(y_mesh+y22).^2))...
        *50.*exp(-k*sqrt((x_mesh-x22).^2 +(y_mesh+y22).^2).*beta);
    v23 = (-gamma2/(2*pi).*(x_mesh-x23)./((x_mesh-x23).^2+(y_mesh-y23).^2))...
        *50.*exp(-k*sqrt((x_mesh-x23).^2 +(y_mesh-y23).^2).*beta);
    vm23 = (gamma2/(2*pi).*(x_mesh+x23)./((x_mesh-x23).^2+(y_mesh+y23).^2))...
        *50.*exp(-k*sqrt((x_mesh-x23).^2 +(y_mesh+y23).^2).*beta);
    v = v1 + vm1 + v21 + vm21 + v22 + vm22 + v23 + vm23;
    
    % Splitting calculation into "2 measurements"
    elevation1 = atan2(y_mesh,x_mesh);
    elevation2 = atan2(y_mesh,(x_mesh-distance));
    v_r1 = u.*cos(elevation1) + v.*sin(elevation1);
    v_r2 = u.*cos(elevation2) + v.*sin(elevation2);
    range1 = sqrt(x_mesh.^2 + y_mesh.^2);
    range2 = sqrt((x_mesh-distance).^2 + y_mesh.^2);
    
    for l=1:size(x_mesh,1)
        for m=1:size(x_mesh,2)
            temp1 = [NaN, range1(l,m), v_r1(l,m), NaN, NaN, elevation1(l,m), NaN, NaN, NaN,...
                NaN, NaN, phi(k)];
            measure_phase_1 = [measure_phase_1; temp1];
            
            temp2 = [NaN, range2(l,m), v_r2(l,m), NaN, NaN, elevation2(l,m), NaN, NaN, NaN,...
                NaN, NaN, phi(k)];
            measure_phase_2 = [measure_phase_2; temp2];
        end
    end
    clear x2* y2* temp* v2* u2* vm2* um2* elev* range*
end

% Implement Doppler intensity into measure_phase_*
measure_phase_1(:,4) = 0.8 + (1.2-0.4).*rand(size(measure_phase_1,1),1);
measure_phase_2(:,4) = 0.8 + (1.2-0.4).*rand(size(measure_phase_2,1),1);

% Put noise on velocity. The velocity values will be manipulated with a
% randomized relative value between a and b
a = 0.9;    b = 1.1;
deviation1 = a + (b-a).*rand(size(measure_phase_1,1),1);
measure_phase_1(:,3) = measure_phase_1(:,3).*deviation1;

deviation2 = a + (b-a).*rand(size(measure_phase_1,1),1);
measure_phase_2(:,3) = measure_phase_2(:,3).*deviation2;

% Change phase from degree to radians
measure_phase_1(:,12) = deg2rad(measure_phase_1(:,12));
measure_phase_2(:,12) = deg2rad(measure_phase_2(:,12));

% Change elevation from radians to degree
measure_phase_1(:,6) = rad2deg(measure_phase_1(:,6));
measure_phase_2(:,6) = rad2deg(measure_phase_2(:,6));
end
