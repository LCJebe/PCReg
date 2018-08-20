clear all
%% READ DATA: model and surface point clouds from saved file
path = 'Data/PointClouds/';
surface_name = 'Surface_DS2_aligned.pcd';
pcSurface = pcread(strcat(path, surface_name));

%% ALIGN: change coordinate system of surface

% specify rotation (in DEG) and translation (in mm)
x_off = 24.5; % 11.7
y_off = 16.7; % 6.1
z_off = 16.2; % 8.3
Rx = 0; % 0
Ry = 0; % 0
Rz = 0; % 0

% build transformation matrix 
R = eul2rotm([Rx, Ry, Rz]/180*pi, 'XYZ');
T = [x_off, y_off, z_off];

A = eye(4);
A(1:3, 1:3) = R;
A(4, 1:3) = T;

% create tform object
tform = affine3d(A);

% transform pointcloud
pcAligned = pctransform(pcSurface,tform);

%% SAVE resulting pointcloud surface as pcd
save_file = strcat(path, surface_name(1:end-4), 'M.pcd');
pcwrite(pcAligned, save_file);