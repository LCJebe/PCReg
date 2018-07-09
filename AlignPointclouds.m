clear all
%% READ DATA: model and surface point clouds from saved file
path = 'Data/PointClouds/';
surface_name = 'Surface_dense.pcd';
pcSurface = pcread(strcat(path, surface_name));
pcModel = pcread(strcat(path, 'GoodCrop.pcd'));

%% ALIGN: change coordinate system of surface

% specify rotation (in DEG) and translation (in mm)
x_off = 10.7; % 11.7
y_off = 6.1; % 6.1
z_off = 8.3; % 8.3
Rx = 70; % 70
Ry = 188; % 188
Rz = 0; % 0

% first: center around the origin
xc = (pcSurface.XLimits(2) + pcSurface.XLimits(1))/2;
yc = (pcSurface.YLimits(2) + pcSurface.YLimits(1))/2;
zc = (pcSurface.ZLimits(2) + pcSurface.ZLimits(1))/2;

A = eye(4);
A(4, 1:3) = -[xc, yc, zc];
tform = affine3d(A);

pcCentered = pctransform(pcSurface,tform);

% build transformation matrix 
R = eul2rotm([Rx, Ry, Rz]/180*pi, 'XYZ');
T = [x_off, y_off, z_off];

A = eye(4);
A(1:3, 1:3) = R;
A(4, 1:3) = T;

% create tform object
tform = affine3d(A);

% transform pointcloud
pcAligned = pctransform(pcCentered,tform);

%% SAVE resulting pointcloud surface as pcd
save_file = strcat(path, surface_name(1:end-4), '_aligned.pcd');
pcwrite(pcAligned, save_file);